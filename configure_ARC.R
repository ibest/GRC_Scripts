

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
    "# numcycles=50",
    "# mapper=bowtie2",
    "# assembler=newbler",
    "# nprocs=40",
    "# format=fastq",
    "# verbose=True",
    "# urt=True",
    "# map_against_reads=False",
    "# assemblytimeout=60",
    "# bowtie2_k=1",
    "# rip=True",
    "# cdna=False",
    "# subsample=1",
    "# maskrepeats=False",
    "Sample_ID\tFileName\tFileType")
  writeLines(config_output,output)
}


write_out_arc_config <- function(reads,sample,output){
  config_output <- c(paste(sample,reads,c("SE","PE1","PE2")[unlist(sapply(c("merged|SE","PE1","PE2"),grep,x=reads))],sep="\t"))
  writeLines(config_output,output)
}


opt <- list(samplesFile="samples.txt",targetsFile="Ref_Sequences/target.txt",readsFolder="02-Cleaned",arcFolder="04-ARC")

### for Rosenblum Whitesands
#opt <- list(samplesFile="S.undulatus_samples.txt",targetsFile="03-targets/S.undulatus/S.undulatus_combined.fasta",readsFolder="02-Cleaned",arcFolder="05-ARC-S.undulatus")
#opt <- list(samplesFile="A.inornata_samples.txt",targetsFile="03-targets/A.inornata/A.inornata_combined.fasta",readsFolder="02-Cleaned",arcFolder="05-ARC-A.inornata")

dir.create(opt$arcFolder)

if(!file.exists(dir(pattern=opt$samplesFile,full.names=TRUE))) stop("Samples file does not exist")
 