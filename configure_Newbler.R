#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("optparse"))
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
  make_option(c("-f", "--file"), type="character", default="samples.txt",
              help="The filename of the sample file [default %default]",
              dest="samplesFile"),
  make_option(c("-c", "--column"), type="character", default="SAMPLE_ID",
              help="Column name from the sample sheet to use as read folder names [default %default]",
              dest="samplesColumn"),  
  make_option(c("-r", "--readFolder"), type="character", default="02-Cleaned",
              help="Directory where the sequence data is stored [default %default]",
              dest="readFolder"),
  make_option(c("-n", "--newblerFolder"), type="character", default="03-NewblerAssemblies",
              help="Directory where to store the newbler results [default %default]",
              dest="newblerFolder"),
  make_option(c("-p", "--processors"), type="integer", default=0,
              help="number of processors to use [defaults to number available]",
              dest="procs"),
  make_option(c("-q", "--newbler_processors"), type="integer", default=10,
              help="number of processors to use in the newbler call [defaults %default]",
              dest="nprocs"),
  make_option(c("-v", "--vector-file"), type="character", default=NULL,
              help="file name with vector sequences in fasta format to provide Newbler [default %default]",
              dest="vector")  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- list(samplesFile="samples.txt",readsFolder="02-Cleaned",newblerFolder="05-Assemblies")
#opt <- list(samplesFile="samples.txt", samplesColumn="SAMPLE_ID", readFolder="04-Screened-extra_human",newblerFolder="05-NewblerAssemblies",  procs=0, nprocs=10)

suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))

######################################################################
## loadSampleFile
## reads in the sample sheet, check for expected format (columns SAMPLE_ID and SEQUENCE_ID)
## then check to make sure sequence reads are available using the specified column
## Parameters
##  file: sample sheet filename
##  reads_folder: path to folder containing reads
##  column: sample sheet column to use that specified folders
"loadSamplesFile" <- function(file, reads_folder,column){
  ## debug
  file = opt$samplesFile; reads_folder = opt$readFolder; column = opt$samplesColumn
  ##
  if ( !file.exists(file) ) {
    write(paste("Sample file",file,"does not exist\n"), stderr())
    stop()
  }  
  ### column SEQUENCE_ID should be the folder name inside of Raw_Folder
  ### column SAMPLE_ID should be the sample name
  ### rows can be commented out with #
  targets <- read.table(file,sep="",header=TRUE,as.is=TRUE)
  if( !all(c("SAMPLE_ID","SEQUENCE_ID") %in% colnames(targets)) ){
    write(paste("Expecting the two columns SAMPLE_ID and SEQUENCE_ID in samples file (tab-delimited)\n"), stderr())
    stop()
  }
  if (any(is.na(match(targets[,column],dir(path=reads_folder))))){
    write(paste(column,"do not match the read data folder structure\n\n"), stderr())
    write(paste(column,"FOUND\n",sep="\t"),stderr())
    write(paste(apply(data.frame(targets[,column],targets[,column] %in% dir(path=reads_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
    stop()
  }
  targets$isDir <- sapply(targets[,column],function(x) file.info(file.path(reads_folder,x))$isdir)
  targets$type <- NA
  for (i in seq.int(to=nrow(targets))){
    if (targets[i,"isDir"]){
      ext <- unique(file_ext(dir(file.path(reads_folder,targets[i,column]),pattern="fastq|sff")))
      if (length(ext) == 0){
        write(paste("Cannot locate fastq or sff file in folder",targets[i,column],"\n"), stderr())
        stop()
      }
      targets$type[i] <- paste(ext,sep="/")
    }
    else {
      ext <- file_ext(grep("fastq|sff",dir(file.path(reads_folder,targets[i,column])),value=TRUE))
      if (length(ext) == 0){
        write(paste(targets[i,column],"is not a fastq or sff file\n"), stderr())
        stop()
      }
      targets$type[i] <- paste(ext,sep="/")
    }
  }
  write(paste("samples sheet contains", nrow(targets), "samples to process",sep=" "),stdout())  
  return(targets)  
}


######################################################################
## prepareCore
##  Set up the numer of processors to use
## 
## Parameters
##  opt_procs: processors given on the option line
##  samples: number of samples
##  targets: number of targets
"prepareCore" <- function(opt_procs){
  # if opt_procs set to 0 then expand to samples by targets
  if( opt_procs == 0 ) opt_procs <- detectCores()
  write(paste("Using",opt_procs,"processors",sep=" "),stdout())
  return(opt_procs)
}

samples <- loadSamplesFile(opt$samplesFile,opt$readFolder,opt$samplesColumn)
procs <- prepareCore(opt$procs)

## create output folder
dir.create(opt$newblerFolder,showWarnings=FALSE,recursive=TRUE)

######################
newblerList <- function(samples, reads_folder, column){
  newbler_list <- list()
  for (i in seq.int(to=nrow(samples))){
    reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="fastq$",full.names=TRUE)
    bt <- lapply(c("_merged|_SE","_PE1|_R1","_PE2|_R2"),grep,x=reads,value=TRUE)
    names(bt) <- c("SE","PE1","PE2")
    bt$sampleFolder=samples[i,column]
    newbler_list[[bt$sampleFolder]] <- bt
  }
  write(paste("Setting up",length(newbler_list),"jobs",sep=" "),stdout())  
  return(newbler_list)
}

newbler <- newblerList(samples,opt$readFolder,opt$samplesColumn)

## run newbler
newbler_out <- mclapply(newbler, function(index){
  dir.create(file.path(opt$newblerFolder,index$sampleFolder))
  try({
    system(paste("runAssembly",
                 "-force -noace -m -sio",
                 "-cpu", opt$nprocs,
                 "-o",file.path(opt$newblerFolder,index$sampleFolder),
                 ifelse(opt$vector == NULL,"",paste("-vt",opt$vector,sep=" ")),
                 paste(index$PE1,collapse=" "),
                 paste(index$PE2,collapse=" "),
                 paste(index$SE,collapse=" "),sep=" "));
  })
},mc.cores=floor(procs/opt$nprocs))

parse_newblerFiles <- function(file){
  lines <- readLines(file)
  # remove extra whitespace
  lines <- gsub("\\s+$", "", lines)
  lines <- lines[!(lines == "")]
  # remove comments
  scomment <- grep("/*",lines,fixed=T)
  ecomment <- grep("*/",lines,fixed=T)
  if (length(scomment) != length(ecomment)){ 
    write(paste("Comment tags aren't matches, file:",file,"\n"), stderr())
    stop()
  }
  lines <- lines[-unlist(apply(cbind(scomment,ecomment),1,function(x) seq.int(x[1],x[2])))]

  procLines <- function(iloop){
    res <- {}
    slist <- which(iloop== "{")
    elist <- which(iloop == "}")
    if (length(slist) > 0 & length(slist) == length(elist)){
      nlist <- split(iloop[unlist(apply(cbind(slist,elist),1,function(x) seq.int(x[1]+1,x[2]-1)))],
                   rep(seq.int(1,length(slist)),times=apply(cbind(slist,elist),1,function(x) (x[2]-x[1]-1))))
      names(nlist) = iloop[slist-1]
      res <- lapply(nlist,function(x) procLines(sub("^\t","",x)))
      iloop = iloop[-unlist(apply(cbind(slist,elist),1,function(x) seq.int(x[1]-1,x[2])))]
    }
    if (length(iloop) > 0){
      iloop <- gsub("^\\s+", "", iloop)
      iloop <- gsub(";|\"| MB", "", iloop)
      ll <- strsplit(iloop,split=" += +")
      nlist <- lapply(ll,"[[",2L)
      names(nlist) <- sapply(ll,"[[",1L)
      res <- c(nlist,res)
    }
    return(res)
  }
  
  # outer list  
  return(procLines(lines))
}

newblertb <- sapply(newbler, function(newb){
  if (!newbler_out[[newb$sampleFolder]]){
    require("Hmisc")
    cfile <- t(data.frame(strsplit(grep("contig",readLines(file.path(opt$newblerFolder,newb$sampleFolder,"454ContigGraph.txt")),value=TRUE),"\t"),stringsAsFactors=F))
    cfile <- data.frame("Contig"=cfile[,2],"Length"=as.numeric(cfile[,3]),"Cov"=as.numeric(cfile[,4]))
    cov <- c(wtd.mean(cfile$Cov,cfile$Length),sqrt(wtd.var(cfile$Cov,cfile$Length)))    
    pfile <- parse_newblerFiles(file.path(opt$newblerFolder,newb$sampleFolder,"454NewblerMetrics.txt"))
    # run data
    areads <- as.numeric(c(pfile$runMetrics$totalNumberOfReads, unlist(strsplit(sub("%","",pfile$consensusResults$readStatus$numAlignedReads),split=" *, *"))))
    abases <- as.numeric(c(pfile$runMetrics$totalNumberOfBases, unlist(strsplit(sub("%","",pfile$consensusResults$readStatus$numAlignedBases),split=" *, *"))))
    rstatus <- as.numeric(c(pfile$consensusResults$readStatus$numberAssembled,
                            pfile$consensusResults$readStatus$numberPartial,
                            pfile$consensusResults$readStatus$numberSingleton,
                            pfile$consensusResults$readStatus$numberRepeat,
                            pfile$consensusResults$readStatus$numberOutlier,
                            pfile$consensusResults$readStatus$numberTooShort))
    passembled <- (sum(rstatus[0:1])/areads[1])*100
    largecontigs <- as.numeric(c(pfile$consensusResults$largeContigMetrics$numberOfContigs,
                                 pfile$consensusResults$largeContigMetrics$numberOfBases,
                                 pfile$consensusResults$largeContigMetrics$avgContigSize,
                                 pfile$consensusResults$largeContigMetrics$N50ContigSize,
                                 pfile$consensusResults$largeContigMetrics$largestContigSize,
                                 unlist(strsplit(sub("%","",pfile$consensusResults$largeContigMetrics$Q40PlusBases),split=" *, *"))))
    allcontigs <- as.numeric(c(pfile$consensusResults$allContigMetrics$numberOfContigs,pfile$consensusResults$allContigMetrics$numberOfBases))
    ndata <- c(areads,abases,rstatus,passembled,largecontigs,allcontigs,cov)
    names(ndata) <- c("totalNumberOfReads","numAlignedReads","numAlignedReadsPercent",
                      "totalNumberOfBases","numAlignedBases","numAlignedReadsBases",
                      "numberAssembled","numberPartial","numberSingleton","numberRepeat","numberOutlier","numberTooShart","assembledPercent",
                      "numLargeContigsAssembled","numLargeBasesAssembled","avgLargeContigSize","N50LargeContigSize","largestContigSize","numQ40PlusBases","Q40PlusBasesPercent",
                      "numAllContigsAssembled","numAllBasesAssembled","meanWeightedCov","sdWeightedCov")
    return(round(ndata,3))
  }
})

newblertb <- t(newblertb)
write.table(newblertb,file.path(opt$newblerFolder,"SummaryNewblerAssemblies.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

