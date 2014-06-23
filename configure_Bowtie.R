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
  make_option(c("-b", "--bowtieFolder"), type="character", default="03-Bowtie",
              help="Directory where to store the bowtie results [default %default]",
              dest="bowtieFolder"),
  make_option(c("-t", "--bowtieTargets"), type="character", default="bowtie_targets.txt",
              help="Path to a bowtie2 build, or a tab delimeted file with name\tbowtie2 index pairs to run bowtie against [default %default]",
              dest="bowtieTarget"),
  make_option(c("-p", "--processors"), type="integer", default=0,
              help="number of processors to use [defaults to number available]",
              dest="procs"),
  make_option(c("-q", "--bowtie_processors"), type="integer", default=10,
              help="number of processors to use in the bowtie call [defaults %default]",
              dest="bprocs"),
  make_option(c("-u", "--extractUnmapped"), action="store_true", default=FALSE,
              help="Extract unmapped reads from the resulting bam file [default %default]",
              dest="extract_unmapped"),
  make_option(c("-m", "--extractMapped"), action="store_true", default=FALSE,
              help="Extract Mapped reads from the resulting bam file [default %default]",
              dest="extract_mapped"),
  make_option(c("-s", "--strict"), action="store_true", default=FALSE,
              help="when extracting mapped reads, use strict (both pairs must map) rules [default %default]",
              dest="strict"),
  make_option(c("-g", "--gzip_extracted"), action="store_true", default=FALSE,
              help="gzip extracted read fastq files [default %default]",
              dest="gzip_extracted"),
  make_option(c("-e", "--extractFolder"), type="character", default="03-Extracted",
              help="if extractUnmapped, and/or extractMapped is TRUE, save resulting fastq to this folder [default %default]",
              dest="screenFolder"),
  make_option(c("-l", "--localmode"), action="store_true", default=FALSE,
              help="use local mode in bowtie2 [default %default]",
              dest="localmode")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- list(samplesFile="samples.txt",readFolder = "02-Cleaned",bowtieFolder="03-Bowtie",bowtieTarget="bowtie_targets.txt",screenFolder="03-Screened",extract_unmapped=TRUE,procs=60)

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
## prepareTargets
##  Prepare the bowtie2 targets to run bowtie2 
## 
## Parameters
##  targets: filename of bowtie2 build, fasta file or text file with multiple targets
"prepareTargets" <- function(targets){
  if (file.exists(paste(targets,"1.bt2",sep="."))){
    ### single target, bowtie2 build exists
    targets_list <- list(c(basename(targets),targets))
  } else if(file_ext(targets) %in% c("fasta","fa","fna")){
    ### single target, need to buld bowtie2 build
    if (!file.exists(paste(sub(".fasta$|.fa$|.fna$","",targets),"1.bt2",sep="."))){
      if (!file.exists(targets)){
        write(paste("Targets file (",targets,") does not exist"))
      }
      write(paste("Preparing bowtie2 indexes for:",targets,"\n"),stdout())
      system(paste("bowtie2-build",targets,sub(".fasta$|.fa$|.fna$","",targets)))      
    }
    targets_list <- list(c(sub(".fasta$|.fa$","",basename(targets)),sub(".fasta$|.fa$|.fna$","",targets)))
  } else if (file.exists(targets)){
    ### multiple targets
    targets_list <- lapply(readLines(targets),function(x) strsplit(x,split="\t")[[1]])
#    if (!all(sapply(targets_list,length) == 2)) {
#      write("Some targets are malformed, this script requires 2 columns (tab separated) per line\n",stderr())
#      stop()
#    }
    for( i in length(targets_list) ) {
      if(file_ext(targets_list[[i]][2]) %in% c("fasta","fa","fna")){
        if (!file.exists(paste(sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2]),"1.bt2",sep="."))){
          write(paste("Preparing bowtie2 indexes for:",targets_list[[i]][2],"\n"),stdout())
          system(paste("bowtie2-build",targets_list[[i]][2],sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2])))
        }
        targets_list[[i]][2] <- sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2])
      }
    }
  } else {
    write(paste("Something wrong with targets file (or table)"),stderr())
    stop()
  }
  write(paste("Found", length(targets_list), "targets to map against",sep=" "),stdout())  
  return(targets_list)
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

######################
"bowtieList" <- function(samples, reads_folder, column, targets){
  bowtie_list <- list()
  for (i in seq.int(to=nrow(samples))){
    reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="fastq$",full.names=TRUE)
    bt <- lapply(c("_merged|_SE","_PE1|_R1","_PE2|_R2"),grep,x=reads,value=TRUE)
    names(bt) <- c("SE","PE1","PE2")
    bt$sampleFolder=samples[i,column]
    for(j in targets){
      bowtie_list[[paste(bt$sampleFolder,j[1],sep="_")]] <- bt
      bowtie_list[[paste(bt$sampleFolder,j[1],sep="_")]]$target_name <- j[1]
      bowtie_list[[paste(bt$sampleFolder,j[1],sep="_")]]$target_path <- j[2]
    }
  }
  write(paste("Setting up",length(bowtie_list),"jobs",sep=" "),stdout())
  return(bowtie_list)
}


samples <- loadSamplesFile(opt$samplesFile,opt$readFolder,opt$samplesColumn)
targets <- prepareTargets(opt$bowtieTarget)
procs <- prepareCore(opt$procs)
bowtie <- bowtieList(samples,opt$readFolder,opt$samplesColumn,targets)

## create output folder
dir.create(opt$bowtieFolder,showWarnings=FALSE,recursive=TRUE)
## run bowtie2
bowtie_out <- mclapply(bowtie, function(index){
  dir.create(file.path(opt$bowtieFolder,index$sampleFolder))
  try({
    system(paste("bowtie2",
        "-I 0 -X 1500",
        ifelse(opt$localmode,"--very-sensitive-local",""),
        "-p", opt$bprocs,
        "--rg-id", index$sampleFolder,
        "--rg", paste("SM",index$sampleFolder,sep=":"),         
        "--rg", paste("PL","illumina",sep=":"),         
        "--rg", paste("LB","whatever",sep=":"),         
        "--rg", paste("PU","whatever",sep=":"),         
        "-x", index$target_path,
        "-1",paste(index$PE1,collapse=","),
        "-2",paste(index$PE2,collapse=","),
        "-U",paste(index$SE,collapse=","),
        #"-S",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"sam",sep=".")),
        "2>",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"out",sep=".")),
        "| samtools view -bS - >", file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),sep=" "));
  })
},mc.cores=floor(procs/opt$bprocs))

## run samtools
samtools_out <- mclapply(bowtie, function(index){
  dir.create(file.path(opt$bowtieFolder,index$sampleFolder))
  try({
    system(paste("samtools sort",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,sep="."))));
    system(paste("samtools index",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),sep=" "));
    system(paste("samtools idxstats",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),">",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"idxstats",sep="."))))
    system(paste("samtools flagstat",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),">",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"flagstat",sep="."))))
  })
},mc.cores=procs)

#####################################################
## write out index tables
targetTables <- sapply(targets,function(tgt){
  filesToRead <- unlist(sapply(file.path(opt$bowtieFolder,unique(samples[,opt$samplesColumn])),dir,pattern=paste(tgt[1],"idxstats",sep="."),full.names=TRUE))
  info <- read.table(filesToRead[1])[,1:2]
  colnames(info) <- c("SequenceID","SequenceLength")
  data <- sapply(filesToRead,function(file){
    tb <- read.table(file)
    values <- rowSums(as.matrix(tb[,3:4]))
    values
  })
  colnames(data) <- basename(colnames(data))
  freq <- round(sweep(data,2,colSums(data),"/")*100,3)
  write.table(cbind(info,data),file.path(opt$bowtieFolder,paste(tgt[1],"summary","reads","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  write.table(cbind(info,freq),file.path(opt$bowtieFolder,paste(tgt[1],"summary","proportions","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  pmapped <- 100-tail(freq,1)
#  as.vector(pmapped)
  names(pmapped) <- colnames(freq)
  round(pmapped,3)
})
targetTables <- data.frame(targetTables)
colnames(targetTables) <- sapply(targets,"[[",1L)

### simple assign by most on target
targetTables <- data.frame(ID=rownames(targetTables),targetTables,assign=colnames(targetTables)[apply(targetTables,1,which.max)])
write.table(targetTables,file.path(opt$bowtieFolder,"SummarySample2Targets.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

#####################################################
## extract unmapped reads
if (opt$extract_unmapped){## Extract Unmapped Reads
  dir.create(file.path(opt$screenFolder))
  extract_out <- mclapply(bowtie, function(index){
    try({
        dir.create(file.path(opt$screenFolder,index$sampleFolder));
        system(paste("samtools view",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_unmapped_reads2.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"unmapped",sep=".")),sep=" "),intern=TRUE);
        #print(paste("samtools view",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_unmapped_reads2.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"unmapped",sep=".")),sep=" "));
    })
  },mc.cores=procs/2)
  extract_out <- strsplit(sapply(extract_out,tail,n=1),split=": |,")
  extract_table <- data.frame(ID=names(bowtie),Records=sapply(extract_out,"[[",2L),PE_pairs=sapply(extract_out,"[[",4L),SE_reads=sapply(extract_out,"[[",6L))
  write.table(extract_table,file.path(opt$screenFolder,"SummaryUnmapped.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

# #####################################################
# ## extract unmapped reads
# if (opt$extract_mapped){## Extract Unmapped Reads
#   dir.create(file.path(opt$screenFolder))
#   extract_out <- mclapply(bowtie, function(index){
#     try({
#       dir.create(file.path(opt$screenFolder,index$sampleFolder));
#       system(paste("samtools view",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_mapped_reads.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"mapped",sep=".")),sep=" "),intern=TRUE);
#       #print(paste("samtools view",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_mapped_reads.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"mapped",sep=".")),sep=" "));
#     })
#   },mc.cores=procs/2)
#   extract_out <- strsplit(sapply(extract_out,tail,n=1),split=": |,")
#   extract_table <- data.frame(ID=names(bowtie),Records=sapply(extract_out,"[[",2L),PE_pairs=sapply(extract_out,"[[",4L),SE_reads=sapply(extract_out,"[[",6L))
#   write.table(extract_table,file.path(opt$screenFolder,"SummaryMapped.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# }

######################################################
## generate VCF files

#java -Xmx6g -jar /mnt/home/msettles/opt/src/GenomeAnalysisTK-2.7-4/GenomeAnalysisTK.jar  -T HaplotypeCaller -R ../../03-targets/A.inornata/A.inornata_combined.fasta -I CP1.A.inornata.bam -o CP1.A.inornata.vcf
# 
# gtk <- "java -Xmx6g -jar /mnt/home/msettles/opt/src/GenomeAnalysisTK-2.7-4/GenomeAnalysisTK.jar"
# if (opt$generate_vcf){## Extract Unmapped Reads
#   vcf_out <- mclapply(bowtie, function(index){
#     try({
#       system(paste(gtk,
#                  "-T HaplotypeCaller",
#                  "-R",paste(index$target_path,"fasta",sep="."),
#                  "-I",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
#                  "-o",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"vcf",sep=".")),
#                  ">",file.path(opt$bowtieFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"gtk","out",sep=".")),sep=" "));    
#     })
#   },mc.cores=procs)
# }
