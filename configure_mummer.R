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
  make_option(c("-r", "--contigFolder"), type="character", default="03-NewblerAssemblies",
              help="Directory where the contig data is stored [default %default]",
              dest="contigFolder"),
  make_option(c('-C', "--contigFileName"), type="character", default="454LargeContigs.fna",
              help="Name of file containing contigs [default %default]", dest="contigFileName"),
  make_option(c("a", "--isArc"),  action="store_true", default=FALSE,
                help="is the folder an ARC folder [default %default]",
                dest="isARC"),
  make_option(c("-m", "--MummerFolder"), type="character", default="04-Mummer",
              help="Directory where to store the mummer results [default %default]",
              dest="mummerFolder"),
  make_option(c("-t", "--mummerTargets"), type="character", default="mummer_targets.txt",
              help="Path to a fasta file, or a tab delimeted file with name\tfasta pairs to run mummer against [default %default]",
              dest="mummerTarget"),
  make_option(c("-p", "--processors"), type="integer", default=0,
              help="number of processors to use [defaults to number available]",
              dest="procs")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- list(samplesFile="samples.txt",readsFolder="02-Cleaned",newblerFolder="05-Assemblies")
#opt <- list(samplesFile="samples.txt", samplesColumn="SAMPLE_ID", readFolder="04-Screened-extra_human",newblerFolder="05-NewblerAssemblies",  procs=0, nprocs=10)

suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("Biostrings"))

######################################################################
## loadSampleFile
## reads in the sample sheet, check for expected format (columns SAMPLE_ID and SEQUENCE_ID)
## then check to make sure sequence reads are available using the specified column
## Parameters
##  file: sample sheet filename
"loadSamplesFile" <- function(file){
  ## debug
#  file = opt$samplesFile; reads_folder = opt$readFolder; column = opt$samplesColumn
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

samples <- loadSamplesFile(opt$samplesFile)
procs <- prepareCore(opt$procs)

## create output folder
dir.create(opt$mummerFolder,showWarnings=FALSE,recursive=TRUE)

"prepareSingleTarget" <- function(targets,path){
  ## targets is a fasta file
  if(file_ext(targets) %in% c("fasta","fa","fna")){
    if (!file.exists(targets)){
      write(paste("Targets file (",targets,") does not exist"))
      stop()
    }
    targets_list <- list(c(sub(".fasta$|.fa$|.fna$","",basename(targets)),targets))
  } else if (file.exists(targets)){
    ### multiple targets
    targets_list <- lapply(readLines(targets),function(x) strsplit(x,split="\t")[[1]])
    if (!all(sapply(targets_list,length) == 2) & !all(sapply(targets_list,length) == 3)) {
      write("Some targets are malformed, this script requires 2 (optionally 3) columns (tab separated) per line\n",stderr())
      stop()
    }
    for( i in length(targets_list) ) {
      if(file_ext(targets_list[[i]][2]) %in% c("fasta","fa","fna")){
        if (!file.exists(targets_list[[i]][2])){
          write(paste("Targets file (",targets_list[[i]][2],") does not exist"))
          stop()
        }
      } else{
        write(paste("Targets file (",targets_list[[i]][2],") is not a fasta, fa or fna file"))
        stop()        
      }
    }
  } else {
    write(paste("Something wrong with targets file (or table)"),stderr())
    stop()
  }
  write(paste("Found", length(targets_list), "targets to map against, copying to mummer folder",sep=" "),stdout()) 
  fas <- sapply(targets_list, function(x) {
      fa <- readDNAStringSet(x[2]);
      names(fa) <- paste(x[1],names(fa),sep="*");
      fa
      })
  writeXStringSet(do.call("c",fas),file.path(path,"targets.fasta"))
  targets_list$combined=file.path(path,"targets.fasta")
  return(targets_list)
}

targets <- prepareSingleTarget(opt$mummerTarget,opt$mummerFolder)
######################
mummerList <- function(samples,contig_folder, contig_file, isARC,targets, column){
  mummer_list <- list()
  if (isARC) samples[,column] <- paste("finished_",samples[,column],sep="")
  for (i in seq.int(to=nrow(samples))){
    contigs <- list(filename=dir(path=file.path(contig_folder,samples[i,column]),pattern=contig_file,full.names=TRUE))
    contigs$sampleFolder=samples[i,column]
    contigs$target = targets
    mummer_list[[contigs$sampleFolder]] <- contigs
  }
  write(paste("Setting up",length(mummer_list),"jobs",sep=" "),stdout())  
  return(mummer_list)
}

if (opt$isARC) opt$contigFileName = "contigs.fasta"
mummer <- mummerList(samples,opt$contigFolder,opt$contigFileName,opt$isARC,targets$combined, opt$samplesColumn)

## run mummer
mummer_out <- mclapply(mummer, function(index){
  dir.create(file.path(opt$mummerFolder,index$sampleFolder))
  try({
    system(paste("mummer",
                 "-L -F -b -mum",
                 index$target,
                 index$filename,
                 ">", file.path(opt$mummerFolder,index$sampleFolder,paste(index$sampleFolder,"mummer",sep=".")),sep=" "),ignore.stderr=TRUE);
  })
},mc.cores=floor(procs))


parse_mummerFiles <- function(file){
  lines <- readLines(file)
  # remove extra whitespace
  lines <- gsub("\\s+$|^\\s+", "", lines)
  lines <- lines[!(lines == "")]

  cindex <- grep ('^>',lines)
  clines <- lines[cindex]  
  clines <- gsub("^> | Len = ","",clines)
  
  clines <- rep(clines,time=diff(c(cindex,length(lines)+1))-1)
  if (length(clines) == 0){
    return(data.frame())
  }
  rev <- grep("Reverse",clines)
  clines <- matrix(unlist(strsplit(sub(" Reverse","",clines),split=" ")),ncol=2,byrow=TRUE)
  colnames(clines) <- c("query_name","query_length")
  lines <- lines[-cindex]
  orient = rep("+",times=length(lines))
  orient[rev] = "-"
  lines <- matrix(unlist(strsplit(lines,split=" +")),ncol=4,byrow=TRUE)
  colnames(lines) <- c("ref_name","pos_ref","pos_query","match_length")
  target <- sub(".","?",sapply(strsplit(lines[,1],split="*",fixed=T),"[[",2L),fixed=T)  
  lines[,1] <- sub(".","?",sapply(strsplit(lines[,1],split="*",fixed=T),"[[",1L),fixed=T)
  clines[,1] <- sub(".","?",clines[,1],fixed=T)
  clines_target <- sapply(strsplit(clines[,1],split="_:_"),"[[",2L)
  
  target_match <- target == clines_target
  
  df <- data.frame(clines,orient,lines,target,clines_target,stringsAsFactors=FALSE)
  df2 <- df[target_match,]
  
  tb <- unlist(
    lapply(split(df2,df2$ref_name),function(x) 
      lapply(split(x,f=x$orient),function(y) 
        sapply(split(y,f=y$query_name),function(z) 
          sum(as.numeric(z[["match_length"]]))/as.numeric(z[["query_length"]][1]))
        )
      )
    )
  tb <- data.frame(matrix(unlist(strsplit(names(tb),split=".",fixed=T)),ncol=3,byrow=T),tb)
  tb <- tb[order(tb[,3],-as.numeric(tb[,4])),]
  tb <- tb[!duplicated(tb[,3]),]
  rownames(tb) <- NULL
  colnames(tb) <- c("REF","ORIENTATION", "QUERY","QUERY COVERAGE")
  tb$"QUERY LENGTH" <- df[match(tb$QUERY,df[,1]),2]
  ### Big dumb rule of 10%
  tb <- tb[tb[["QUERY COVERAGE"]] >= 0.1,]
  tb$BREAK = NA
  df.tmp <- df2[match(paste(tb$QUERY,tb$ORIENTATION),paste(df2$query_name,df2$orient)),]
  tb$BREAK[df.tmp$pos_ref == 1 & df.tmp$pos_query != 1] <- df.tmp$pos_query [df.tmp$pos_ref == 1 & df.tmp$pos_query != 1]

  #df <- df[(paste(df$query_name,df$orient) %in% paste(tb$QUERY,tb$ORIENTATION)),c(1:7)]
  df <- df[(paste(df$query_name) %in% paste(tb$QUERY)),c(1:7)]
  df$query_name <- sub("?",".",df$query_name,fixed=T)
  df$ref_name <- sub("?",".",df$ref_name,fixed=T)

  tb$REF <- sub("?",".",tb$REF,fixed=T)
  tb$QUERY <- sub("?",".",tb$QUERY,fixed=T)
  list(tb=tb,df=df)
  
}


mummertb <- sapply(mummer, function(index)
  if (!mummer_out[[index$sampleFolder]]){
#    cat(index$sampleFolder,"\n")
    pmummer <- parse_mummerFiles(file.path(opt$mummerFolder,index$sampleFolder,paste(index$sampleFolder,"mummer",sep=".")))
    df <- pmummer$df
    write.table(df,file.path(opt$mummerFolder,index$sampleFolder,paste(index$sampleFolder,"mummer.txt",sep=".")),sep="\t",row.names=F,col.names=T,quote=F)
    tb <- pmummer$tb
    if (length(tb) > 0){
      write.table(tb,file.path(opt$mummerFolder,index$sampleFolder,paste(index$sampleFolder,"target.txt",sep=".")),sep="\t",row.names=F,col.names=T,quote=F)
      fa <- readDNAStringSet(index$filename)
      nms <- sapply(strsplit(names(fa),split=" "),"[[",1L)
      tb <- tb[match(nms,tb$QUERY),]
      names(fa) <- paste(tb$REF,index$sampleFolder,names(fa),sep="-")
      targets$combined= NULL
      if (all(sapply(targets,length) == 3)){
        keep <- as.logical(sapply(targets,"[[",3L)[match(tb$REF,sapply(targets,"[[",1L))])
        keep[is.na(keep)] <- FALSE
        fa <- fa[keep]
        tb <- tb[keep,]
      }
      for (i in seq.int(1,length(fa))){
        if (!is.na(tb[i,"BREAK"])){
          fa2 <- subseq(fa[i],as.numeric(tb[i,"BREAK"]))
          fa[i] <- subseq(fa[i],1, as.numeric(tb[i,"BREAK"])-1)
          nms <- sapply(strsplit(names(fa[i]),split=" "),"[[",1L)
          names(fa)[i] <- paste(nms,".1 SPLIT_CONTIG REGION 1",sep="")
          names(fa2) <- paste(nms,".2 SPLIT_CONTIG REGION 2",sep="")
          fa <- c(fa,fa2)
        }
      }
      fa <- fa[order(names(fa))]
      writeXStringSet(fa,file.path(opt$mummerFolder,index$sampleFolder,paste(index$sampleFolder,"screened.fasta",sep=".")))    
    }
})

