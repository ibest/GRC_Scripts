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
