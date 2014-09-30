
opt <- list(
    samplesFile="samples.txt",
    targetsFile="transposon.fa",
    readsFolder="02-Cleaned",
    arcFolder="03-ARC")


arc_header <- function(ref,output){
  config_output <- c(
    "## Name=value pairs:",
    "## reference: contains reference sequences in fasta format",
    "## numcycles: maximum number of times to try remapping",
    "## mapper: the mapper to use (blat/bowtie2)",
    "## assembler: the assembler to use (newbler/spades)",
    "## nprocs: number of cores to use",
    "## format: fastq or fastq, all must be the same",
    "##",
    "## Columns:",
    "## Sample_ID:Sample_ID",
    "## FileName: path for fasta/fastq file",
    "## FileType: PE1, PE2, or SE",
    "## FileFormat: fasta or fastq",
    paste("# reference=",ref,sep=""),
    "# numcycles=20",
    "# mapper=bowtie2",
    "# assembler=newbler",
    "# nprocs=60",
    "# format=fastq",
    "# verbose=True",
    "# urt=True",
    "# map_against_reads=False",
    "# assemblytimeout=60",
    "# bowtie2_k=1",
    "# rip=True",
    "# cdna=False",
    "# subsample=0.1",
    "# maskrepeats=False",
    "Sample_ID\tFileName\tFileType")
  writeLines(config_output,output)
}


write_out_arc_config <- function(reads,sample,output){
  config_output <- c(paste(sample,reads,c("SE","PE1","PE2")[unlist(sapply(c("merged|SE","PE1","PE2"),grep,x=reads))],sep="\t"))
  writeLines(config_output,output)
}


### for Rosenblum Whitesands
#opt <- list(samplesFile="S.undulatus_samples.txt",targetsFile="03-targets/S.undulatus/S.undulatus_combined.fasta",readsFolder="02-Cleaned",arcFolder="05-ARC-S.undulatus")
#opt <- list(samplesFile="A.inornata_samples.txt",targetsFile="03-targets/A.inornata/A.inornata_combined.fasta",readsFolder="02-Cleaned",arcFolder="05-ARC-A.inornata")

dir.create(opt$arcFolder)

if(!file.exists(dir(pattern=opt$samplesFile,full.names=TRUE))) stop("Samples file does not exist")

samples <- read.table(opt$samplesFile,sep="",header=T,as.is=T)
if(!file.exists(opt$targetsFile)) stop("Targets file does not exist")

zz <- file(file.path(opt$arcFolder,"ARC_config.txt"), "w")  # open an output file connection
arc_header(normalizePath(opt$targetsFile),zz)

for (i in unique(samples$SAMPLE_ID)){
  read2arc <- dir(file.path(getwd(),opt$readsFolder,i),pattern="fastq",full.names=TRUE)  
  write_out_arc_config(read2arc,i,zz)         
}

close(zz)

