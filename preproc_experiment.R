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
  make_option(c("-d", "--directory"), type="character", default="00-RawData",
              help="Directory where the raw sequence data is stored [default %default]",
              dest="Raw_Folder"),
  make_option(c("-q", "--quality"), type="integer", default=24,
              help="Quality score to use during lucy trimming [default %default]",
              dest="qual"),
  make_option(c("-m", "--miniumumLength"), type="integer", default=150,
              help="Discard reads less then minimum length [default %default]",
              dest="minL"),
  make_option(c("-o", "--overlap"), type="integer", default=700,
              help="Overlap parameter for flash [default %default]",
              dest="overlap"),
  make_option(c("-O", "--skip-overlap"), action="store_true", default=FALSE,
              help="do not perform the overlapping using flash [default %default]",
              dest="noOverlap"),
  make_option(c("-p", "--processors"), type="integer", default=1,
              help="number of processors to use [default %default]",
              dest="procs"),
  make_option(c("-s", "--skip-duduplicates"), action="store_true", default=FALSE,
              help="do not perform the deduplication step [default %default] ",
              dest="skip_dedup"),
  make_option(c("-c", "--contaminants-folder"), type="character", default=NULL,
              help="folder name with contaminant sequences in fasta format [default %default]",
              dest="contaminants"),
  make_option(c("-v", "--vector-folder"), type="character", default=NULL,
              help="folder name with vector sequences in fasta format [default %default]",
              dest="vector"),
  make_option(c("-a", "--polyA"), action="store_true", default=FALSE,
              help="perform polyA trimming in seqyclean [default %default]",
              dest="polyA"),
  make_option(c("--i64"), action="store_true",default=FALSE,
              help="input read Q scores are offset by 64 [default %default]",
              dest="i64")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ShortRead"))
suppressPackageStartupMessages(library("parallel"))
####################################################
### FUNCTIONS
####################################################
screen_duplicates <- function(r,o,d,s){
  paste("screen_duplicates_PE.py",ifelse(s,"-s",""),"-d", r, "-o", o, ">>", file.path(d,"preprocessing_output.txt"), sep=" ")
#  paste("screen_duplicates_PE_sra.py","-d", r, "-o", o, ">>", file.path(d,"preprocessing_output.txt"), sep=" ")
}

seqyclean_illumina <- function(r1,r2,o,minL=150, q=24,polyA,folder, sample, i64) {
  i64_param = ""
  if (i64){
   i64_param="-i64"
  }
  vc_param = ""
  if(file.exists(file.path(folder,"contaminants.fa"))){
    vc_param=paste(vc_param,"-c",file.path(getwd(),folder,"contaminants.fa"),sep=" ")
  }
  if(file.exists(file.path(folder,"vector.fa"))){
    vc_param=paste(vc_param,"-v",file.path(getwd(),folder,"vector.fa"),sep=" ")
  }
	paste("seqyclean --ow -qual", q, q, i64_param, vc_param,ifelse(polyA,"-polyat",""),"-minimum_read_length",minL,"--new2old_illumina -1",r1,"-2",r2,"-o",o, ">>", file.path(folder,sample,"preprocessing_output.txt"),sep=" ")
}

join_reads <- function(r1,r2,o,overlap=275,d){
#	paste("flash --allow-outies --max-overlap=",overlap," --output-prefix=",o," ",r1," ",r2, " >> ", file.path(d,"preprocessing_output.txt"),sep="")
    paste("flash --max-overlap=",overlap," --output-prefix=",o," ",r1," ",r2, " >> ", file.path(d,"preprocessing_output.txt"),sep="")
}

link_illumina <- function(se1,se2,r1,r2,o){
  require("ShortRead")
  output = ""
  ## first merge the 2 SE files
  if (!is.na(se1) & !is.na(se2)){
      se <- file.path(dirname(se1),"merged_SE_files.fastq")
      system(paste("cat",se1,se2,">",se,sep=" "))
#     if (!file.exists(se)){
#         fq <- readFastq(c(se1,se2))
#         writeFastq(fq,se)
# 	}
  } else{
    se = se1
  }
  output <- paste("mv -f",file.path(se),paste(o,"merged_SE.fastq",sep="_"),";","mv -f",file.path(r1),paste(o,"notcombined_PE1.fastq",sep="_"),";","mv -f",file.path(r2),paste(o,"notcombined_PE2.fastq",sep="_"),";",sep=" ")
  output
}


seqyclean_454 <- function(sff,o,minL=225,q=24,polyA,folder,sample){
    vc_param = ""
  if(file.exists(file.path(folder,"contaminants.fa"))){
    vc_param=paste(vc_param,"-c",file.path(getwd(),folder,"contaminants.fa"),sep=" ")
  }
  if(file.exists(file.path(folder,"vector.fa"))){
    vc_param=paste(vc_param,"-v",file.path(getwd(),folder,"vector.fa"),sep=" ")
  }
	paste("seqyclean -qual",q,q,vc_param,ifelse(polyA,"-polyat",""),"-minimum_read_length",minL,"-454",sff,"-o",o, ">>", file.path(folder,sample,"preprocessing_output_454.txt"),sep=" ")
}
link_454 <- function(sff,o){
	paste("mv -f",sff,paste(o,"sff",sep="."),sep=" ")
}

final_report_fun <- function(f,o){
  paste("read_info.py", "-d",f,">>",file.path(o,"preprocessing_output_final_report.txt"))
}

get_phiX <- function(){
  genbank_id <- "NC_001422"
  URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                 genbank_id, "&rettype=fasta&retmode=text", 
                 sep = "")
  res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
  phix <- DNAStringSet(paste(res[-1],collapse=""))
  names(phix) <- "contaminant_PhiX"
  return(phix)
}

#qual <- opt$qual
#minL <- 120
#overlap <- opt$overlap
process_sample <- function(folder,sample,Raw_Folder,Clean_Folder,Final_Folder,qual,polyA,minL,overlap,noOverlap,i64,skipd){
  write(paste(sample,":Processing folder ",folder,sep=""),stdout())
  if(file.info(file.path(Raw_Folder,folder))$isdir){ ## ILLUMINA FOLDER, EXPECT PAIRED READS
    
    if(!file.exists(file.path(Clean_Folder, sample))) dir.create(file.path(Clean_Folder, sample),recursive=TRUE,showWarnings=FALSE)
    
    output <- file.path(Clean_Folder,sample,paste(sample,sep="_"))
    if (skipd){
        write(paste(sample,":\tmerging any reads in RawData",sep=""),stdout())
    }else{
        write(paste(sample,":\tde-duplicating reads",sep=""),stdout())        
    }

    system(screen_duplicates(file.path(Raw_Folder,folder),output,file.path(Clean_Folder,sample),skipd))

    ## second use seqyclean to remove contaminant, adapters and trim for quality
    Read1 <- paste(output,"nodup_PE1.fastq.gz",sep="_")
    Read2 <- paste(output,"nodup_PE2.fastq.gz",sep="_")
    output <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),sep="_"))
    write(paste(sample,":\trunning seqyclean",sep=""),stdout())

    seqyclean_cmd <- seqyclean_illumina(Read1, Read2, output, minL=minL, q=qual, polyA,Clean_Folder, sample, i64)
    system(seqyclean_cmd)

    ## third use flash to join overlapping paired-end reads
    Read1 <- paste(output,"_PE1.fastq",sep="")
    Read2 <- paste(output,"_PE2.fastq",sep="")
    output <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),sep="_"))
    write(paste(sample,":\tjoining reads",sep=""),stdout())

    SE1 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),"SE.fastq",sep="_"))
    SE2 <- NA
    if (!noOverlap){
        system(join_reads(Read1, Read2,output,overlap=overlap,file.path(Clean_Folder,sample)))
        ## link the final files to another folder
        SE2 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".extendedFrags.fastq",sep=""),sep="_"))
    
        Read1 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".notCombined_1.fastq",sep=""),sep="_"))
        Read2 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".notCombined_2.fastq",sep=""),sep="_"))
    }

    output <- file.path(getwd(),Final_Folder,sample,sample)
    
    if(!file.exists(file.path(Final_Folder, sample))) dir.create(file.path(Final_Folder,sample),recursive=TRUE,showWarnings=FALSE)
    write(paste(sample,":\tMoving final files to ",file.path(Final_Folder,sample),sep=""),stdout())
    system(link_illumina(SE1,SE2,Read1,Read2,output))
    write(paste(sample,":\tFinished",sep=""),stdout())

    ## end preprocessing
  } else { ## 454 READS
    qual454 <- 20
    minL454 <- 250
    if(!file.exists(file.path(Clean_Folder, sample))) dir.create(file.path(Clean_Folder,sample),recursive=TRUE,showWarnings=FALSE)
    SFF <- file.path(Raw_Folder,folder)
    folder <- sub(".sff","",folder)
    output <- file.path(getwd(),Clean_Folder,sample,paste(sample,paste("q",qual454,"min",minL454,sep=""),sep="_"))
    seqyclean_cmd <- seqyclean_454(SFF,output,minL=minL454,q=qual454,polyA,Clean_Folder, sample)    
    write(paste(sample,":\trunning 454 seqyclean",sep=""),stdout())
    system(seqyclean_cmd)
    SFF <- paste(output,"sff",sep=".")
    output <- file.path(getwd(),Final_Folder,sample,sample)
    if(!file.exists(file.path(Final_Folder, sample))) dir.create(file.path(Final_Folder,sample),recursive=TRUE,showWarnings=FALSE)
    write(paste(sample,":\tcreating 454 links to final files in ",file.path(Final_Folder,sample),sep=""),stdout())
    system(link_454(SFF,output))
    write(paste(sample,":\t454 Finished",sep=""),stdout())
  }	
}


##########################################

##########################################
## SETUP PROCESSING
## test
# opt <- list(samplesFile="samples.txt",Raw_Folder="00-RawData",qual=24,minL=150,overlap=275,procs=1)

if ( !file.exists(opt$samplesFile) ) {
  write(paste("Sample file",opt$samplesFile,"does not exist\n"), stderr())
  stop()
}

### opt$samplesFile$SEQUENCE_ID should be the folder name inside of Raw_Folder
### opt$samplesFile$SAMPLE_ID should be the sample name
targets <- read.table(opt$samplesFile,sep="\t",header=T,as.is=T)

if (detectCores() < opt$procs){
  write(paste("number of cores specified (",opt$procs,") is greater than the number of cores available (",detectCores(),")",sep=" "),stdout())
  stop()
}

##########################################
## directory structure
## assumes Illumina data is under a folder, named opt$samplesFile$SEQUENCE_ID and 454 sff files are just files in the folder
Raw_Folder <- opt$Raw_Folder
## where preprocessing clean results are saved
Clean_Folder <- "01-Clean_Merge"
if(!file.exists(file.path(Clean_Folder))) dir.create(file.path(Clean_Folder),recursive=TRUE,showWarnings=FALSE)

## final preprocess file are linked (from Clean_Folder) to this folder
Final_Folder <- "02-Cleaned"
if(!file.exists(file.path(Final_Folder))) dir.create(file.path(Final_Folder),recursive=TRUE,showWarnings=FALSE)

if (any(is.na(match(targets$SEQUENCE_ID,dir(path=Raw_Folder))))){
  write("samples file does not match the raw data folder structure\nExpecting two columns (tab delimited) with column headings SEQUENCE_ID and SAMPLE_ID", stderr())
  stop()
}

##########################################
### main function

if (!is.null(opt$contaminants)){
  if(!file.info(opt$contaminants)$isdir){
    write(paste("Parameter contaminant must be a directory",sep=" "), stderr())
    stop()
  }
  files <- dir(opt$contaminants,pattern="fa|fasta",full.names=TRUE)
  write(paste("found",length(files),"contaminant files to join",sep=" "),stdout())
  contaminants <- readDNAStringSet(files)
  contaminants <- c(contaminants,get_phiX()) ## add in phiX
  writeXStringSet(contaminants,file.path(Clean_Folder,"contaminants.fa"))
} else {
  contaminants <- get_phiX() ## add in phiX
  writeXStringSet(contaminants,file.path(Clean_Folder,"contaminants.fa"))
}

if (!is.null(opt$vector)){
  if(!file.info(opt$vector)$isdir){
    write(paste("Parameter vector must be a directory",sep=" "), stderr())
    stop()
  }
  files <- dir(opt$vector,pattern="fa|fasta",full.names=TRUE)
  write(paste("found",length(files),"vector files to join",sep=" "),stdout())
  vector <- readDNAStringSet(files)
  writeXStringSet(vector,file.path(Clean_Folder,"vector.fa"))
}

write(paste("samples sheet contains", nrow(targets), "samples to process",sep=" "),stdout())


mclapply(seq.int(1,nrow(targets)), function(index){
  folder <- targets$SEQUENCE_ID[index]
  sample <- targets$SAMPLE_ID[index]
  try({process_sample(folder,sample,Raw_Folder,Clean_Folder,Final_Folder,opt$qual,opt$polyA,opt$minL,opt$overlap,opt$noOverlap,opt$i64,opt$skip_dedup)})
  write(paste(sample,":\tcreating report of final files in ",file.path(Final_Folder,sample),sep=""),stdout())
  try({system(final_report_fun(file.path(Final_Folder,sample),file.path(Clean_Folder,sample)))})
},mc.cores=opt$procs)


##########################################
### Generate Read Report
write(paste("Generating Final Preprocessing Report for all samples",sep=" "),stdout())

system(paste("preproc_report -f", opt$samplesFile,"-c",Clean_Folder,"-d",Final_Folder,sep=" "))

write("Finished processing samples",stdout())

