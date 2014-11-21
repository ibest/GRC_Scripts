###### R functions involved in Preprocessing

## run the python application screen_duplicates
# r dirctory that contains the reads
# o output directory
# d output results
screen_duplicates <- function(r,o,d){
  paste("screen_duplicates_PE.py","-d", r, "-o", o, ">>", file.path(d,"preprocessing_output.txt"), sep=" ")
  #  paste("screen_duplicates_PE_sra.py","-d", r, "-o", o, ">>", file.path(d,"preprocessing_output.txt"), sep=" ")
}

## run seqyclean on illumina data
# r1, path to read1
# r2, path to read2
# o, output location
# minL, minimum length
# q, quality trim param
# folder, folder containing contaminant and vector fasta files (contaminants.fa and/or vector.fa)
# sample, sample name
# i64, is data in illumina 64 format
seqyclean_illumina <- function(r1,r2,o,minL=150, q=24,folder, sample, i64) {
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
  paste("seqyclean --ow -qual", q, q, i64_param, vc_param,"-minimum_read_length",minL,"--new2old_illumina -1",r1,"-2",r2,"-o",o, ">>", file.path(folder,sample,"preprocessing_output.txt"),sep=" ")
}

## run the application flash to join overlapping paired end reads
join_reads <- function(r1,r2,o,overlap=275,d){
  paste("flash --max-overlap=",overlap," --output-prefix=",o," ",r1," ",r2, " >> ", file.path(d,"preprocessing_output.txt"),sep="")
}

## generate final linked files in cleaned folder with report
link_illumina <- function(se1,se2,r1,r2,o){
  require("ShortRead")
  ## first merge the 2 SE files
  se <- file.path(dirname(se1),"merged_SE_files.fastq")
  if (!file.exists(se)){
    fq <- readFastq(c(se1,se2))
    writeFastq(fq,se)
  }
  paste("ln -sf",file.path("../..",se),paste(o,"merged_SE.fastq",sep="_"),";ln -sf",file.path("../..",r1),paste(o,"notcombined_PE1.fastq",sep="_"),";ln -sf",file.path("../..",r2),paste(o,"notcombined_PE2.fastq",sep="_"),";",sep=" ")
}

## run seqyclean on 454 data
# path to sff file
# o, output location
# minL, minimum length
# q, quality trim param
# folder, folder containing contaminant and vector fasta files (contaminants.fa and/or vector.fa)
# sample, sample name
seqyclean_454 <- function(sff,o,minL=225,q=24,folder,sample){
  vc_param = ""
  if(file.exists(file.path(folder,"contaminants.fa"))){
    vc_param=paste(vc_param,"-c",file.path(getwd(),folder,"contaminants.fa"),sep=" ")
  }
  if(file.exists(file.path(folder,"vector.fa"))){
    vc_param=paste(vc_param,"-v",file.path(getwd(),folder,"vector.fa"),sep=" ")
  }
  paste("seqyclean -qual",q,q,vc_param,"-minimum_read_length",minL,"-454",sff,"-o",o, ">>", file.path(folder,sample,"preprocessing_output_454.txt"),sep=" ")
}
link_454 <- function(sff,o){
  paste("ln -sf",sff,paste(o,"sff",sep="."),sep=" ")
}

final_report_fun <- function(f,o){
  paste("read_info.py", "-d",f,">>",file.path(o,"preprocessing_output_final_report.txt"))
}

## download phiX from NCBI for screen Illumina data, always perform
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

##### Actual Preprocessing Workflow
process_sample <- function(folder,sample,Raw_Folder,Clean_Folder,Final_Folder,qual,minL,overlap,i64){
  write(paste(sample,":Processing folder ",folder,sep=""),stdout())
  if(file.info(file.path(Raw_Folder,folder))$isdir){ ## ILLUMINA FOLDER, EXPECT PAIRED READS
    
    if(!file.exists(file.path(Clean_Folder, sample))) dir.create(file.path(Clean_Folder, sample),recursive=TRUE,showWarnings=FALSE)
    
    output <- file.path(Clean_Folder,sample,paste(sample,sep="_"))
    write(paste(sample,":\tde-duplicating reads",sep=""),stdout())
    
    system(screen_duplicates(file.path(Raw_Folder,folder),output,file.path(Clean_Folder,sample)))
    
    ## second use seqyclean to remove contaminant, adapters and trim for quality
    Read1 <- paste(output,"nodup_PE1.fastq.gz",sep="_")
    Read2 <- paste(output,"nodup_PE2.fastq.gz",sep="_")
    output <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),sep="_"))
    write(paste(sample,":\trunning seqyclean",sep=""),stdout())
    
    seqyclean_cmd <- seqyclean_illumina(Read1, Read2, output, minL=minL, q=qual, Clean_Folder, sample, i64)
    system(seqyclean_cmd)
    
    ## third use flash to join overlapping paired-end reads
    Read1 <- paste(output,"_PE1.fastq",sep="")
    Read2 <- paste(output,"_PE2.fastq",sep="")
    output <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),sep="_"))
    write(paste(sample,":\tjoining reads",sep=""),stdout())
    
    system(join_reads(Read1, Read2,output,overlap=overlap,file.path(Clean_Folder,sample)))
    
    ## link the final files to another folder
    SE1 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,sep=""),"SE.fastq",sep="_"))
    SE2 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".extendedFrags.fastq",sep=""),sep="_"))
    
    Read1 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".notCombined_1.fastq",sep=""),sep="_"))
    Read2 <- file.path(Clean_Folder,sample,paste(sample,"nodup",paste("q",qual,"min",minL,".notCombined_2.fastq",sep=""),sep="_"))
    output <- file.path(getwd(),Final_Folder,sample,sample)
    
    if(!file.exists(file.path(Final_Folder, sample))) dir.create(file.path(Final_Folder,sample),recursive=TRUE,showWarnings=FALSE)
    write(paste(sample,":\tcreating links to final files in ",file.path(Final_Folder,sample),sep=""),stdout())
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
    seqyclean_cmd <- seqyclean_454(SFF,output,minL=minL454,q=qual454,Clean_Folder, sample)    
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

