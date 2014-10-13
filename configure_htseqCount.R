#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("parallel"))

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
    make_option(c("-b", "--mappingFolder"), type="character", default="03-BWA",
                help="Directory where to store the mapping results [default %default]",
                dest="mappingFolder"),
    make_option(c("-t", "--mappingTargets"), type="character", default="mapping_targets.txt",
                help="Path to a gtf file, or tab delimeted file with [target name]\t[target fasta]\t[target gtf] to run mapping against [default %default]",
                dest="mappingTarget"),
    make_option(c("-h", "--htseqFolder"), type="character", default="04-HTseqCounts",
                help="Directory where to store the HtSeq-cont results [default %default]",
                dest="htseqFolder"),
    make_option(c("-r", "--order"), metavar = "POS", type="character", default="name",
                help="'pos' or 'name'. Sorting order of <alignment_file> [default: %default]. Paired-end sequencing data must be sorted either by position or by read name, and the sorting order must be specified. Ignored for single-end data.",
                dest="order"),
    make_option(c("-s", "--stranded"), metavar = "STRANDED", type="character", default="yes",
                help="whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' [default: %default]. 'reverse' means 'yes' with reversed strand interpretation",
                dest="stranded"),
    make_option(c("-a", "--minaqual"), metavar = "MINAQUAL", type="integer", default=10,
                help="skip all reads with alignment quality lower than the given minimum value [default: %default]",
                dest="minaqual"),
    make_option(c("-y", "--type"), metavar = "FEATURETYPE", type="character", default="exon",
                help="feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)",
                dest="type"),
    make_option(c("-i", "--idattr"), metavar = "IDATTR", type="character", default="gene_id",
                help="GFF attribute to be used as feature ID (default, suitable for Ensembl GTF files: %default)",
                dest="idattr"),
    make_option(c("-m", "--mode"), metavar = "MODE", type="character", default="union",
                help="mode to handle reads overlapping more than one feature (choices: union, intersection-strict, intersection-nonempty; default: %default)",
                dest="mode"),                
    make_option(c("-p", "--processors"), type="integer", default=0,
                help="number of processors to use [defaults to number available]",
                dest="procs")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

######################################################################
## prepareCore
##    Set up the numer of processors to use
## 
## Parameters
##	opt_procs: processors given on the option line
##	samples: number of samples
##	targets: number of targets
"prepareCore" <- function(opt_procs){
    # if opt_procs set to 0 then expand to samples by targets
    if( opt_procs == 0 ) opt_procs <- detectCores()
    write(paste("Using up to",opt_procs,"processors",sep=" "),stdout())
    return(opt_procs)
}

procs <- prepareCore(opt$procs)

######################################################################
## loadSampleFile
## reads in the sample sheet, check for expected format ('SAMPLE_ID' or column parameter)
## Parameters
##  file: sample sheet filename, column name for the sample ids
"loadSamplesFile" <- function(file,column){
    ##
    if ( !file.exists(file) ) {
        write(paste("Sample file",file,"does not exist\n"), stderr())
        stop("Quiting")
    }  
    ### column SEQUENCE_ID should be the folder name inside of Raw_Folder
    ### column SAMPLE_ID should be the sample name
    ### rows can be commented out with #
    targets <- read.table(file,sep="",header=TRUE,as.is=TRUE)
    if( !(column %in% colnames(targets)) ){
        write(paste("Expecting", column,"in the header of samples file\n",sep=" "), stderr())
        stop("Quiting")
    }
    write(paste("Samples sheet contains", nrow(targets), "samples to process",sep=" "),stdout())  
    return(targets)  
}

samples <- loadSamplesFile(opt$samplesFile,opt$samplesColumn)

######################################################################
## loadDataFile
## provided the sample sheet,
## then check to make sure sequence reads are available using the specified column
## Parameters
##	file: sample sheet filename
##	data_folder: path to folder containing reads
##	column: sample sheet column to use that specified folders
##  raw: if processing raw data, folderids use a differnt column identifier 'raw'
"loadDataFile" <- function(samples, data_folder,column,raw=NA){
    ##
    if (!is.na(raw)){  ### RAW READS IN RAW READ FOLDER, USES INFO FROM A DIFFERENT COLUMN IN THE SAMPLESHEET
        if (any(is.na(samples[,raw] %in% dir(path=data_folder)))){
            write(paste(raw,"do not match the raw data folder structure\n\n"), stderr())
            write(paste(raw,"FOUND\n",sep="\t"),stderr())
            write(paste(apply(data.frame(samples[,raw],samples[,raw] %in% dir(path=data_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
            write("\n",stderr())
            stop("Quitting")
        }
        samples$isDir <- sapply(samples[,raw],function(x) file.info(file.path(data_folder,x))$isdir)
        samples$type <- NA
        for (i in seq.int(to=nrow(samples))){
            if (samples[i,"isDir"]){
                ext <- unique(file_ext(dir(file.path(data_folder,samples[i,raw]),pattern="fastq|sff")))
                if (length(ext) == 0){
                    write(paste("Cannot locate fastq or sff file in folder",samples[i,raw],"\n"), stderr())
                    stop()
                }
                samples$type[i] <- paste(ext,sep="/")
            }
            else {
                ext <- file_ext(grep("fastq|sff",dir(file.path(data_folder,samples[i,raw])),value=TRUE))
                if (length(ext) == 0){
                    write(paste(samples[i,raw],"is not a fastq or sff file\n"), stderr())
                    stop()
                }
                samples$type[i] <- paste(ext,sep="/")
            }
        }
        write(paste("data found, ready to process",sep=" "),stdout())    
        return(samples)	
    }
    #### NOT raw data, processed data
    if (any(is.na(samples[,column] %in% dir(path=data_folder)))){
        write(paste(column,"do not match the data folder structure\n\n"), stderr())
        write(paste(column,"FOUND\n",sep="\t"),stderr())
        write(paste(apply(data.frame(samples[,column],samples[,column] %in% dir(path=data_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
        write("\n",stderr())        
        stop()
    }
    # remove any duplicate rownames
    samples <- samples[!duplicated(samples[,column]),]
    samples$isDir <- sapply(samples[,column],function(x) file.info(file.path(data_folder,x))$isdir)
    samples$type <- NA
    write(paste("data found, ready to process",sep=" "),stdout())    
    return(samples)    
}

targets <- loadDataFile(samples, opt$mappingFolder,opt$samplesColumn,raw=NA)

######################################################################
## prepareTargets
##    Prepare the mapping targets 
## 
## Parameters
##	targets: filename of targets builds, fasta file or text file with multiple targets
"prepareTargets" <- function(targets, algorithm=NA){
    ### single target, indexes exist
    if ((algorithm == "bowtie" & file.exists(paste(targets,"rev.2.bt2",sep="."))) | (algorithm == "bwa" & file.exists(paste(targets,"bwt",sep=".")))){
        ### single target, bowtie2 build exists
        targets_list <- list(c(basename(targets),targets))
    } else if( file_ext(targets) %in% c("fasta","fa","fna") ){
        ### single target, need to build indexes
        if (!file.exists(targets)){
            write(paste("Targets file (",targets,") does not exist"), stderr())
            stop("Quiting")
        }
        if (algorithm == "bowtie"){
            if(!file.exists(paste(sub(".fasta$|.fa$|.fna$","",targets),"rev.2.bt2",sep="."))){
                write(paste("Preparing bowtie2 indexes for:",targets,"\n"),stdout())
                res <- system(paste("bowtie2-build",targets,sub(".fasta$|.fa$|.fna$","",targets)),ignore.stdout=T, ignore.stderr=T)
                if (res != 0){
                    write(paste("Failed building Bowtie2 indexes for (",targets,") "), stderr())
                    stop("Quiting")
                }
            }
            targets_list <- list(c(sub(".fasta$|.fa$|.fna$","",basename(targets)),sub(".fasta$|.fa$|.fna$","",targets)))	
        } else if (algorithm == "bwa"){
            if(!file.exists(paste(targets,"bwt",sep="."))){
                write(paste("Preparing bwa indexes for:",targets,"\n"),stdout())
                res <- system(paste("bwa index",targets),ignore.stdout=T, ignore.stderr=T) 
                if (res != 0){
                    write(paste("Failed building BWA indexes for (",targets,") "), stderr())
                    stop("Quiting")
                }		
            }
            targets_list <- list(c(basename(targets),targets))					
        }
    } else if( file_ext(targets) %in% c("gtf","gff") ){
        ### single target, gff files
        if (!file.exists(targets)){
            write(paste("gtf [or gff] file (",targets,") does not exist"), stderr())
            stop("Quiting")
        }
        targets_list <- list(c(basename(targets),"NA",targets))
    } else if (file.exists(targets)){
        ### multiple targets
        targets_list <- lapply(readLines(targets),function(x) strsplit(x,split="\t")[[1]])
        #	 Assume first column is name, second is the fasta file, remaining columns are ignored
        for( i in seq.int(length(targets_list)) ) {
            if(file_ext(targets_list[[i]][2]) %in% c("fasta","fa","fna")){
                if (!file.exists(targets_list[[i]][2])){
                    write(paste("Targets file (",targets_list[[i]][2],") does not exist"), stderr())
                    stop("Quiting")
                }
                if (algorithm == "bowtie"){
                    if(!file.exists(paste(sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2]),"rev.2.bt2",sep="."))){
                        write(paste("Preparing bowtie2 indexes for:",targets_list[[i]][2],"\n"),stdout())
                        res <- system(paste("bowtie2-build",targets_list[[i]][2],sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2])),ignore.stdout=T, ignore.stderr=T)
                        if (res != 0){
                            write(paste("Failed building Bowtie2 indexes for (",targets_list[[i]][2],") "), stderr())
                            stop("Quiting")
                        }
                    }
                    targets_list[[i]][2] <- sub(".fasta$|.fa$|.fna$","",targets_list[[i]][2])
                } else if (algorithm == "bwa"){
                    if(!file.exists(paste(targets_list[[i]][2],"bwt",sep="."))){
                        write(paste("Preparing bwa indexes for:",targets_list[[i]][2],"\n"),stdout())
                        res <- system(paste("bwa index",targets_list[[i]][2]),ignore.stdout=T, ignore.stderr=T) 
                        if (res != 0){
                            write(paste("Failed building BWA indexes for (",targets_list[[i]][2],") "), stderr())
                            stop("Quiting")
                        }		
                    }
                    targets_list[[i]][2] <- targets_list[[i]][2]			 
                }
            } else if(file_ext(targets_list[[i]][3]) %in% c("gtf","gff")){
                if (!file.exists(targets_list[[i]][3])){
                    write(paste("Targets file (",targets_list[[i]][3],") does not exist"), stderr())
                    stop("Quiting")
                }
                targets_list[[i]][2] <- targets_list[[i]][3]
            } else {
                write(paste("Something wrong with targets file (or table)"),stderr())
                stop("Quiting")
            }
        }
    } else {
        write(paste("Something wrong with targets file (or table)"),stderr())
        stop("Quiting")
    }
    write(paste("Found", length(targets_list), "targets",sep=" "),stdout())	
    return(targets_list)
}

targets <- prepareTargets(opt$mappingTarget)


######################
"processingList" <- function(samples, data_folder, column, targets, type="sequence"){
    mapping_list <- list()
    for (i in seq.int(to=nrow(samples))){
        if (type == "sequence"){
            files <- dir(path=file.path(data_folder,samples[i,column]),pattern="fastq$",full.names=TRUE)
            map <- lapply(c("TEST","_merged|_SE","_PE1|_R1","_PE2|_R2"),grep,x=files,value=TRUE)
            names(map) <- c("TEST","SE","PE1","PE2")
            map$sampleFolder=samples[i,column]
            for(j in targets){
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]] <- map
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_name <- j[1]
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_path <- j[2]
            }
        } else if (type == "bam"){
            if (length(targets) > 1){
                for (j in targets){
                    files <- dir(path=file.path(data_folder,samples[i,column]),pattern=paste(j[1],".+bam$",sep=""),full.names=TRUE)
                    map <- list("readsorted"= grep("byreadid",files,value=TRUE),
                                "positionsorted"= grep("byreadid",files,invert=TRUE,value=TRUE))
                    map$sampleFolder=samples[i,column]
                    ## Right now only works for 1 target
                    mapping_list[[paste(map$sampleFolder,j[1],sep="_")]] <- map
                    mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_name <- j[1]
                    mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_path <- j[3]                    
                }
                write(paste("Sorry can only process 1 target at the most, please provide the script with a gtf file"),stderr())
                stop("Quiting")
            } else {
                files <- dir(path=file.path(data_folder,samples[i,column]),pattern="bam$|sam$",full.names=TRUE)
                map <- list("readsorted"= grep("byreadid",files,value=TRUE),
                            "positionsorted"= grep("byreadid",files,invert=TRUE,value=TRUE))
                map$sampleFolder=samples[i,column]
                ## Right now only works for 1 target
                j = targets[[1]]
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]] <- map
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_name <- j[1]
                mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_path <- j[3]
            }
        }
    }
    write(paste("Setting up",length(mapping_list),"jobs",sep=" "),stdout())
    return(mapping_list)
}

processing <- processingList(samples,opt$mappingFolder,opt$samplesColumn,targets, type="bam")

## create output folder
dir.create(opt$htseqFolder,showWarnings=FALSE,recursive=TRUE)

htseq_out <- mclapply(processing, function(index){
    dir.create(file.path(opt$htseqFolder,index$sampleFolder),showWarnings=FALSE)
    try({
        call = paste("htseq-count -f bam",
                     "-r", opt$order,
                     "-s", opt$stranded,
                     "-a", opt$minaqual,
                     "-t", opt$type,
                     "-i", opt$idattr,
                     "-m", opt$mode,
                     ifelse(opt$order=="pos",index$positionsorted,index$readsorted),
                     index$target_path,
                     "2>",file.path(opt$htseqFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"out",sep=".")),
                     ">", file.path(opt$htseqFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"counts",sep=".")),sep=" ");
        system(call)
    })
},mc.cores=procs)

if (!all(sapply(htseq_out, "==", 0L))){
    write(paste("Something went wrong with htseq-count, some jobs failed"),stderr())
    stop("Quiting")
}

#####################################################
## write out summary tables
htseqTables <- sapply(targets,function(tgt){
    print(paste("Generating output for target:",tgt[1]))
    filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$htseqFolder,x,paste(x,tgt[1],"counts",sep="."))))
    #	filesToRead <- unlist(sapply(file.path(opt$mappingFolder,unique(samples[,opt$samplesColumn])),dir,pattern=paste(tgt[1],"idxstats",sep="."),full.names=TRUE))
    info <- lapply(filesToRead,read.table,sep="\t",as.is=TRUE)
    names <- info[[1]][,1]
    statidx <- grep("__",names)
    stat = sapply(info,function(x) x[statidx,2])
    info = sapply(info,function(x) x[-statidx,2])
    
    htseq_data <- data.frame("Reads in feature"=colSums(info),"Reads NOT in feature"=stat[1,],"Reads ambiguous"=stat[2,],"Reads too low qual"
                    =stat[3,],"Percent Assigned To Feature"=colSums(info)/(colSums(info)+colSums(stat)),"Number of Features"=nrow(info),"Number of 0 count features"=apply(info,2,function(x)sum(x == 0)))
    write.table(htseq_data,file.path(opt$htseqFolder,paste(tgt[1],"summary","txt",sep=".")),row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
#    htseq_data
})
