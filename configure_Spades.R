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
  make_option(c("-n", "--spadesFolder"), type="character", default="03-SpadesAssemblies",
              help="Directory where to store the spades results [default %default]",
              dest="spadesFolder"),
  make_option(c("-p", "--processors"), type="integer", default=0,
              help="number of processors to use [defaults to number available]",
              dest="procs"),
  make_option(c("-q", "--spades_processors"), type="integer", default=10,
              help="number of processors to use in the spades call [defaults %default]",
              dest="nprocs"),
  make_option(c("-e", "--error-correction"), action="store_true", default=FALSE,
              help="Perform spades error correction prior to assembly [default %default]",
              dest="errorCorrect"),
  make_option(c("-l", "--large_contig_size"), type="integer", default=500,
              help="size of contig considered 'large' [defaults %default]",
              dest="largeContig")
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
dir.create(opt$spadesFolder,showWarnings=FALSE,recursive=TRUE)

######################
spadesList <- function(samples, reads_folder, column){
  spades_list <- list()
  for (i in seq.int(to=nrow(samples))){
    reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="fastq$",full.names=TRUE)
    bt <- lapply(c("_merged|_SE","_PE1|_R1","_PE2|_R2"),grep,x=reads,value=TRUE)
    names(bt) <- c("SE","PE1","PE2")
    bt$sampleFolder=samples[i,column]
    spades_list[[bt$sampleFolder]] <- bt
  }
  write(paste("Setting up",length(spades_list),"jobs",sep=" "),stdout())  
  return(spades_list)
}

spades <- spadesList(samples,opt$readFolder,opt$samplesColumn)

if (length(spades$PE1) > 1 | length(spades$SE > 1)) stop("Sorry but configure_Spades.R does not handle more than one fastq read set")
## run newbler
spades_out <- mclapply(spades, function(index){
  if(file.exists(file.path(opt$spadesFolder,index$sampleFolder))) suppressWarnings(unlink(file.path(opt$spadesFolder,index$sampleFolder),recursive=TRUE))
  dir.create(file.path(opt$spadesFolder,index$sampleFolder))
  try({
    call <- paste("spades.py",
                 "--careful",
#                 "-k 21,33,55,77,99,127",
                 "-t", opt$nprocs,
                 ifelse(opt$errorCorrect,"","--only-assembler"), ## default --only-assembler
                 "-o",file.path(opt$spadesFolder,index$sampleFolder),
                 ifelse(length(index$PE1),paste(
                    "-1",index$PE1[1],"-2",index$PE2[1],sep=" "),""),
                 ifelse(length(index$SE),paste(
                    "-s",index$SE[1],sep=" "),""),
                 ">", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"spades","out",sep=".")),sep=" ");
    res <- system(call);        
    if (res == 0 & file.exists(file.path(opt$spadesFolder,index$sampleFolder,"scaffolds.fasta"))){ # OK
        library(Biostrings)
        contig_seqs <- readDNAStringSet(file.path(opt$spadesFolder,index$sampleFolder,"contigs.fasta"))
        writeXStringSet(contig_seqs[width(contig_seqs) > opt$largeContig],file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"LargeContigs.fna",sep=".")))
        writeXStringSet(contig_seqs,file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"AllContigs.fna",sep=".")))

        scaffold_seqs <- readDNAStringSet(file.path(opt$spadesFolder,index$sampleFolder,"scaffolds.fasta"))
        writeXStringSet(scaffold_seqs[width(scaffold_seqs) > opt$largeContig],file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds.fna",sep=".")))
        scaffolds <- file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds.fna",sep="."))        
        res <- system(paste("bwa index",scaffolds),ignore.stdout=T, ignore.stderr=T)
        if (res == 0){ ## BWA Map reads back to scaffolds
            res = res2 = res3 = 0        
            if(length(index$PE1)) try({
                call <- paste("bwa mem",
                                  "-M", # Mark shorter split hits as secondary (for Picard compatibility).
#                                  "-L 0,0", # penalty for 5' or 3' end clipping 
                                  "-t", opt$nprocs,
                                  "-R", paste("'@RG",
                                              paste("ID",index$sampleFolder,sep=":"),         
                                              paste("SM",index$sampleFolder,sep=":"),    			 
                                              paste("PL","ILLUMINA",sep=":"),				 
                                              paste("LB","whatever",sep=":"),				 
                                              paste("PU","whatever",sep=":"),
                                              paste("DS","Paired",sep=":"), "'",sep="\t"),			 
                                  scaffolds,
                                  paste(index$PE1,collapse=","),
                                  paste(index$PE2,collapse=","),
                                  "2>",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","out",sep=".")),
                                  "| samtools view -bS -F 0x100 - 2> /dev/null >", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")),sep=" ");
                res<-system(call);
                system(paste("samtools view  -H", 
                             file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")),
                             "| head -n -1 > ",
                             file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","header",sep=".")),
                             "2> /dev/null",sep=" "))
            })
            if(length(index$SE)) try({
                res2<-system(paste("bwa mem",
                                   "-M", # Mark shorter split hits as secondary (for Picard compatibility).
                                   "-t", opt$nprocs,
                                   "-R", paste("'@RG",
                                               paste("ID",index$sampleFolder,sep=":"),         
                                               paste("SM",index$sampleFolder,sep=":"),    			 
                                               paste("PL","ILLUMINA",sep=":"),				 
                                               paste("LB","whatever",sep=":"),				 
                                               paste("PU","whatever",sep=":"),
                                               paste("DS","Paired",sep=":"),"'", sep="\t"),			 
                                   scaffolds,
                                   paste(index$SE,collapse=","),
                                   "2>",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","out",sep=".")),
                                   "| samtools view -bS -F 0x100 - 2> /dev/null >", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep=".")),sep=" "));
                system(paste("samtools view  -H", 
                             file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")),
                             "| head -n -1 > ",
                             file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","header",sep=".")),
                             "2> /dev/null",sep=" "))
            })
            ## reheader
            if (file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep="."))))
                file.remove(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")))
            if (length(index$PE1) & length(index$SE)){
                res3 <- system(paste("samtools cat -h", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","header",sep=".")), "-o",
                                     file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),
                                     ifelse(file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep="."))),
                                            file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")),""),
                                     ifelse(file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep="."))),
                                            file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep=".")),""),
                                     "2> /dev/null",sep=" "))
            } else if (file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")))){
                res3 <- system(paste("samtools reheader", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","header",sep=".")),
                                     file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")),
                                     ">", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),
                                     "2> /dev/null",sep=" "))
            } else if (file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep=".")))){
                res3 <- system(paste("samtools reheader", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","header",sep=".")),
                                     file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep=".")),
                                     ">", file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),
                                     "2> /dev/null",sep=" "))
            } else {
                res3 = 1
            }
            if (file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep="."))))
                file.remove(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","PE","bam",sep=".")))     
            if (file.exists(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep="."))))
                file.remove(file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","SE","bam",sep=".")))     
        }   
    }
    try({
        res4 <- system(paste("samtools sort",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds",sep=".")),"2> /dev/null",sep=" "));
        res4 <- res4 & system(paste("samtools index",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),"2> /dev/null",sep=" "));
        res4 <- res4 & system(paste("samtools idxstats",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),">",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","idxstats",sep=".")),"2> /dev/null",sep=" "))
        res4 <- res4 & system(paste("samtools flagstat",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","bam",sep=".")),">",file.path(opt$spadesFolder,index$sampleFolder,paste(index$sampleFolder,"Scaffolds","flagstat",sep=".")),"2> /dev/null",sep=" "))
    })    
  })
  return(as.integer(res | res2 | res3 | res4))
},mc.cores=floor(procs/opt$nprocs))


#####################################################
## write out report tables

calcN50 <- function(contigLengths){
    contigs<-rev(sort(contigLengths))
    contigs[cumsum(contigs) >= sum(contigs)/2][1]    
}

calcCov <- function(coverages,contigLengths){
    require("Hmisc")
    c("meanWeightedCov" = wtd.mean(coverages,contigLengths),"sdWeightedCov" = sqrt(wtd.var(coverages,contigLengths)))    
}

filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$spadesFolder,x,paste(x,"Scaffolds","flagstat",sep="."))))
data.flagstat <- sapply(filesToRead,function(file){
    values <- readLines(file)
    values <- sapply(strsplit(values,split=" + 0",fixed=T),"[[",1L)
    as.numeric(values)
})
rownames(data.flagstat) <- c("totalNumberOfReads","secondary","supplementary","duplicates","numMappedReads","ReadsPaired","read1","read2","ProperlyPaired","itselfandmate","Singletons","mappedAcrossContigs","mapChrQ5")
data.flagstat = rbind(data.flagstat[c("totalNumberOfReads","numMappedReads"),],"numMappedReadsPercent"=data.flagstat["numMappedReads",]/data.flagstat["totalNumberOfReads",],data.flagstat[c("ReadsPaired","ProperlyPaired"),],"ProperlyPairedPercent"=data.flagstat["ProperlyPaired",]/data.flagstat["ReadsPaired",],data.flagstat[c("Singletons","mappedAcrossContigs"),])


filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$spadesFolder,x,paste(x,"AllContigs","fna",sep="."))))
data.AllContigs <- sapply(filesToRead,function(file){
    IDs <- sub("^>","",grep(pattern = "^>",readLines(file),value = TRUE))
    values <- strsplit(IDs,split="_length_|_cov_|_ID_")
    names = sapply(values,"[[",1L)
    lengths = as.numeric(sapply(values,"[[",2L))
    cov = as.numeric(sapply(values,"[[",3L))
    c("numAllContigsAssembled"=length(names),"numAllBasesAssembled"=sum(lengths),"avgAllContigSize"=mean(lengths),"N50AllContigSize"=calcN50(lengths))
})

filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$spadesFolder,x,paste(x,"LargeContigs","fna",sep="."))))
data.LargeContigs <- sapply(filesToRead,function(file){
    IDs <- sub("^>","",grep(pattern = "^>",readLines(file),value = TRUE))
    values <- strsplit(IDs,split="_length_|_cov_|_ID_")
    names = sapply(values,"[[",1L)
    lengths = as.numeric(sapply(values,"[[",2L))
    cov = as.numeric(sapply(values,"[[",3L))
    c("numLargeContigsAssembled"=length(names),"numLargeBasesAssembled"=sum(lengths),"avgLargeContigSize"=mean(lengths),"N50LargeContigSize"=calcN50(lengths))
})


filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$spadesFolder,x,paste(x,"Scaffolds","idxstats",sep="."))))
data.idxstats <- sapply(filesToRead,function(file){
    idxstats <- head(read.table(file,as.is=TRUE),-1)
    values <- strsplit(idxstats$V1,split="_length_|_cov_|_ID_")
    names <- sapply(values,"[[",1L)
    lengths <- as.numeric(sapply(values,"[[",2L))
    cov <- as.numeric(sapply(values,"[[",3L))
    idxstats <- data.frame("Scaffold_FULL"=idxstats$V1,"Scaffold_ID"=names,"Scaffold_Length"=lengths,"Scaffold_Cov"=cov,"Scaffold_MappedReads"=idxstats$V3)
    write.table(idxstats,paste(file,"txt",sep="."),sep="\t",col.names=T,row.names=F,quote=F)
    c("numScaffoldContigsAssembled"=length(names),"numScaffoldBasesAssembled"=sum(lengths),"avgScaffoldContigSize"=mean(lengths),"N50ScaffoldContigSize"=calcN50(lengths),"largestContigSize"=max(lengths),calcCov(cov,lengths))
})

finalSummary <- t(rbind(data.flagstat,data.AllContigs,data.LargeContigs,data.idxstats))

#################################################################################

write.table(finalSummary,file.path(opt$spadesFolder,"SummarySpadesAssemblies.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

