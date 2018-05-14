

### Spencer data, link over fastq files

Raw_Folder <- "00-RawData"

data_folder <- "/mnt/home/uirig/user_data/Spencer_T/Ovine-22-Transcriptomes-Sept-2012/RawData"

sub_folder <- c("Berkeley-Sept-2012","Berkeley-Oct-2012")

sf <- sub_folder[2]

dir(file.path(data_folder,sf),pattern="Sample")

files2Link <- dir(dir(file.path(data_folder,sf,"Project_IBEST_GRC"),pattern="Sample",full.names=TRUE),pattern=".fastq.gz",full.names=T)
s_folder <- sapply(strsplit(files2Link,split="/"),"[[",11L)
s_names <- paste(sf,sapply(strsplit(files2Link,split="/"),"[[",12L),sep="_")

file.symlink(files2Link,file.path(Raw_Folder,s_folder,s_names))


## Bifido data, generate assemblies

samples <- read.table("samples.txt",sep="",header=T,as.is=T)

seq_file <- "02-Cleaned"
assm_file <- "03-Newbler_Assembles"

if(file.exists(assm_file)) unlink(assm_file,recursive=TRUE)
dir.create(assm_file,recursive=TRUE,showWarnings=FALSE)

library(parallel)
runNewbler <- function(sample,cov=NA,cpu=5){
  output_dir <- file.path(assm_file,sample)
  reads <- paste(dir(file.path(seq_file,sample),pattern="fastq",full.names=TRUE),collapse=" ")
  if (is.na(cov)){
    paste("runAssembly -force","-cpu", cpu,"-o",output_dir,reads,sep=" ")
  }else{
    paste("runAssembly -force","-cpu", cpu,"-e",cov,"-o",output_dir,reads,sep=" ")    
  }
}


############# RUN NEWBLER
mclapply(samples$SAMPLE_ID,function(x) system(runNewbler(x)),mc.cores=60)

############# GET COVERAGE, RUN NEWBLER AGAIN
mclapply(samples$SAMPLE_ID,function(x){
  tb <- read.table(file.path(assm_file,x,"454ContigGraph.txt"),sep="\t",as.is=T)
  tb <- tb[substring(tb$V2,1,6) == "contig",]
  cov <- weighted.mean(as.numeric(tb$V4),as.numeric(tb$V3))
  system(runNewbler(x,cov=cov))
  writeLines(runNewbler(x,cov=cov),con=file.path(assm_file,x,"NewblerCall.txt"))
},mc.cores=60)


########### MUMMER 
#########################################
library(parallel)
assm_file <- "03-Newbler_Assembles"

samples <- read.table("samples.txt",sep="",header=T,as.is=T)

patric_db <- "/mnt/home/uirig/user_data/genomes/patric"
patric_species <- dir(patric_db,pattern="Bifidobacterium|Enterococcus")

nucmer_pairs <- cbind(rep(samples$SAMPLE_ID,each = length(patric_species)),rep(patric_species, times=length(samples$SAMPLE_ID)))
nucmer_pairs <- apply(nucmer_pairs,1,as.list)

#ref_Bifido <- "/mnt/home/uirig/user_data/genomes/patric/Bifidobacterium_longum/Bifidobacterium_longum.fna"
#ref_Entero <- "/mnt/home/uirig/user_data/genomes/patric/Enterococcus_faecalis/Enterococcus_faecalis.fna"

map_file <- "04-Nucmer"
if(file.exists(map_file)) unlink(map_file,recursive=TRUE)
dir.create(map_file,recursive=TRUE,showWarnings=FALSE)

runMummer <- function(sample,ref){
  if(!file.exists(file.path(map_file,sample))) dir.create(file.path(map_file,sample))
  output_dir <- file.path(map_file,sample,paste(ref,"mummer.txt",sep="."))
  ref_file <- file.path("/mnt/home/uirig/user_data/genomes/patric",ref,paste(ref,"fna",sep="."))
  if (file.exists(ref_file))
    paste("mummer","-b -F -c -L", ref_file,file.path(assm_file,sample,"454AllContigs.fna"),">>",output_dir) 
}

runNucmer <- function(sample,ref){
  if(!file.exists(file.path(map_file,sample))) dir.create(file.path(map_file,sample))
  
  output_dir <- file.path(map_file,sample,ref)
  ref_file <- file.path("/mnt/home/uirig/user_data/genomes/patric",ref,paste(ref,"fna",sep="."))
  if (file.exists(ref_file))
    paste("nucmer","-p",output_dir, ref_file ,file.path(assm_file,sample,"454AllContigs.fna")) 
}

tmp <- mclapply(nucmer_pairs,function(x){
  system(runNucmer(x[[1]],x[[2]]))
  system(runMummer(x[[1]],x[[2]]))
},mc.cores=60)



library(parallel)
library(Biostrings)

assm_file <- "03-Newbler_Assembles"
samples <- read.table("samples.txt",sep="",header=T,as.is=T)

for (sample in samples$SAMPLE_ID){
  fa <- readDNAStringSet(file.path(assm_file,sample,"454LargeContigs.fna"))
  nms <- sapply(strsplit(names(fa),split=" +"),"[[",1L)
  info <- read.table(file.path(assm_file,sample,"454ContigGraph.txt"),sep="\t",as.is=TRUE)
  info <- info[match(nms,info$V2),]
  print(table(round(as.numeric(info[,4]),-1)))

}

######################## BRICE RESULTING CONTIGS
library(Biostrings)
path <- "/scratch/shunter/tamias"
path <- "/mnt/home/shunter/Tamias_Exome"
contigs <- DNAStringSet(do.call(c,unname(sapply(sapply(dir(path=path,pattern="finished",full.names=T),dir,pattern="processedcontigs.fa",full.names=T),function(x) as.character(readDNAStringSet(x))))))
nms <- names(contigs)
#dt <- as.data.frame(do.call(rbind,strsplit(nms,split=" +|=|(_:_)"))[,c(1,2,4,6,8)],stringsAsFactors=FALSE)
#dt$V4 <- as.numeric(dt$V4)
#dt$V5 <- as.numeric(dt$V5)
#colnames(dt) <- c("sample","target","contig","length","numreads")
dt <- as.data.frame(do.call(rbind,lapply(strsplit(nms,split="\t|(_:_)"),function(x) x[c(1,2,5)])),stringsAsFactors=FALSE)
colnames(dt) <- c("sample","target","type")
keep <- which(dt$type %in% c("overlapped","single"))
dt <- dt[keep,]
contigs <- contigs[keep,]
#targets <- readDNAStringSet("/scratch/msettles/Sarver/Tamias_Exome/target.fasta")
targets <- readDNAStringSet("/mnt/home/uirig/GRC_projects/Sarver/Tamias_Exome/2013-08-23-Tamias_revised_targets/revised_merged_targets.fasta")
tnms <- sapply(strsplit(names(targets),split="_:_"),"[[",2L)
target_size <- rowSums(alphabetFrequency(targets)[,c("A","C","T","G")])
names(target_size) <- names(targets)

#unique_targets <- names(table(tnms))
tb <- sapply(tapply(dt$target,dt$sample,table),function(x) x[match(unique_targets,names(x))])
rownames(tb) <- unique_targets

writeXStringSet(contigs,"AllResults.fa")

source("http://webpages.uidaho.edu/msettles/Rcode/blat.R")
Alltargets <- readPSL("All2Target.psl")
Alltargets$ident <- pslIdent(Alltargets)
Alltargets$cov <- pslCoverage(Alltargets)

ord <- order(Alltargets$qName,Alltargets$tName,Alltargets$tStart)
Alltargets <- Alltargets[ord,]

ContigQuery <- sapply(strsplit(Alltargets$qName,split="_:_"),"[[",2L)
ContigTarget <- sapply(strsplit(Alltargets$tName,split="_:_"),"[[",2L)
identID <- paste(Alltargets$qName,Alltargets$qStart,Alltargets$qEnd,Alltargets$matches)
which(identID %in% (identID[duplicated(identID)])


ContigQueryNew <- sapply(strsplit(AlltargetNewer$qName,split="_:_"),"[[",2L)
ContigTargetNew <- sapply(strsplit(AlltargetNewer$tName,split="_:_"),"[[",2L)

qnamesplit <- split()


ContigQueryFinal <- sapply(strsplit(FinalTargets$qName,split="_:_"),"[[",2L)
ContigTargetFinal <- sapply(strsplit(FinalTargets$tName,split="_:_"),"[[",2L)
      