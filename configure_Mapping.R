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
	make_option(c("-r", "--readFolder"), type="character", default="02-Cleaned",
							help="Directory where the sequence data is stored [default %default]",
							dest="readFolder"),
	make_option(c("-M", "--mappingAlgorithm"), type="character", default="bwa",
							help="Mapping algorithm to use, supported types are 'bowtie' and 'bwa' [default %default]",
							dest="mappingAlgorithm"),
	make_option(c("-l", "--localmode"), action="store_true", default=FALSE,
							help="use local mode in bowtie2 [default %default]",
							dest="localmode"),	
	make_option(c("-b", "--mappingFolder"), type="character", default=NA,
							help="Directory where to store the mapping results [default '03-[mappingAlgorithm]']",
							dest="mappingFolder"),
	make_option(c("-t", "--mappingTargets"), type="character", default="mapping_targets.txt",
	                        help="Path to a fasta file, or tab delimeted file with [target name]\t[target fasta]\t[target gtf, optional] to run mapping against [default %default]",
            	            dest="mappingTarget"),
	make_option(c("-p", "--processors"), type="integer", default=0,
							help="number of processors to use [defaults to number available]",
							dest="procs"),
	make_option(c("-q", "--mapping_processors"), type="integer", default=10,
							help="number of processors to use in the mapping call [defaults %default]",
							dest="mprocs"),
	make_option(c("-n", "--sortByReadID"), action="store_true", default=FALSE,
	                        help="When sorting bam files, sort by read ID (samtools -n option), for compatability with htseq-count [default %default]",
	                        dest="sortByReadID"),
	make_option(c("-i", "--ignoreSingles"), action="store_true", default=FALSE,
	                        help="Ignore any single-end files, for compatability with htseq-count [default %default]",
	                        dest="ignoreSingles"),
	make_option(c("-u", "--extractUnmapped"), action="store_true", default=FALSE,
							help="Extract unmapped reads from the resulting bam file [default %default]",
							dest="extract_unmapped"),
	make_option(c("-m", "--extractMapped"), action="store_true", default=FALSE,
							help="Extract Mapped reads from the resulting bam file [default %default]",
							dest="extract_mapped"),
	make_option(c("-e", "--extractFolder"), type="character", default="04-Extracted",
							help="if extractUnmapped, and/or extractMapped is TRUE, save resulting fastq to this folder [default %default]",
							dest="screenFolder"),	
	make_option(c("-s", "--strict"), action="store_true", default=FALSE,
							help="when extracting mapped reads, use strict (both pairs must map) rules [default %default]",
							dest="strict"),
	make_option(c("-g", "--gzip_extracted"), action="store_true", default=FALSE,
							help="gzip extracted read fastq files [default %default]",
							dest="gzip_extracted")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

if (!(opt$mappingAlgorithm %in% c("bowtie","bwa")))
    stop("Mapping algorithm parameter must be one of 'bowtie' or 'bwa'")

if (is.na(opt$mappingFolder)){
    if (opt$mappingAlgorithm == "bowtie")
        opt$mappingFolder = "03-Bowtie"
    else if (opt$mappingAlgorithm == "bwa")
        opt$mappingFolder = "03-BWA"
    else
        stop("Error in setting mapping folder")
}
######################################################################
## loadSampleFile
## reads in the sample sheet, check for expected format (columns SAMPLE_ID and SEQUENCE_ID)
## then check to make sure sequence reads are available using the specified column
## Parameters
##	file: sample sheet filename
##	reads_folder: path to folder containing reads
##	column: sample sheet column to use that specified folders
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
##	Prepare the mapping targets 
## 
## Parameters
##	targets: filename of targets builds, fasta file or text file with multiple targets
"prepareTargets" <- function(targets, algorithm){
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
			}
		}
	} else {
		write(paste("Something wrong with targets file (or table)"),stderr())
		stop("Quiting")
	}
	write(paste("Found", length(targets_list), "targets to map against",sep=" "),stdout())	
	return(targets_list)
}

######################################################################
## prepareCore
##	Set up the numer of processors to use
## 
## Parameters
##	opt_procs: processors given on the option line
##	samples: number of samples
##	targets: number of targets
"prepareCore" <- function(opt_procs){
	# if opt_procs set to 0 then expand to samples by targets
	if( opt_procs == 0 ) opt_procs <- detectCores()
	write(paste("Using",opt_procs,"processors",sep=" "),stdout())
	return(opt_procs)
}

######################
"mappingList" <- function(samples, reads_folder, column, targets, algorithm){
	mapping_list <- list()
	for (i in seq.int(to=nrow(samples))){
		reads <- dir(path=file.path(reads_folder,samples[i,column]),pattern="fastq$",full.names=TRUE)
		map <- lapply(c("TEST","_merged|_SE","_PE1|_R1","_PE2|_R2"),grep,x=reads,value=TRUE)
		names(map) <- c("TEST","SE","PE1","PE2")
        if (opt$ignoreSingles) map$SE=character(0)
		map$sampleFolder=samples[i,column]
		for(j in targets){
			mapping_list[[paste(map$sampleFolder,j[1],sep="_")]] <- map
			mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_name <- j[1]
			mapping_list[[paste(map$sampleFolder,j[1],sep="_")]]$target_path <- j[2]
		}
	}
	write(paste("Setting up",length(mapping_list),"jobs",sep=" "),stdout())
	return(mapping_list)
}


samples <- loadSamplesFile(opt$samplesFile,opt$readFolder,opt$samplesColumn)
targets <- prepareTargets(opt$mappingTarget, opt$mappingAlgorithm)
procs <- prepareCore(opt$procs)
mapping <- mappingList(samples,opt$readFolder,opt$samplesColumn,targets,opt$mappingAlgorithm)

## create output folder
dir.create(opt$mappingFolder,showWarnings=FALSE,recursive=TRUE)

if (opt$mappingAlgorithm == "bowtie"){
    ## run bowtie2
    bowtie_out <- mclapply(mapping, function(index){
    	dir.create(file.path(opt$mappingFolder,index$sampleFolder),showWarnings=FALSE)
    	try({
    		system(paste("bowtie2",
    				"-I 0 -X 1500",
    				ifelse(opt$localmode,"--very-sensitive-local",""),
    				"-p", opt$mprocs,
    				"--rg-id", index$sampleFolder,
    				"--rg", paste("SM",index$sampleFolder,sep=":"),				 
    				"--rg", paste("PL","ILLUMINA",sep=":"),				 
    				"--rg", paste("LB","whatever",sep=":"),				 
    				"--rg", paste("PU","whatever",sep=":"),				 
    				"-x", index$target_path,
                    ifelse(length(index$PE1),paste(
        				"-1",paste(index$PE1,collapse=","),
        				"-2",paste(index$PE2,collapse=","),
                        sep=" "),""),
                    ifelse(length(index$SE),paste(
        				"-U",paste(index$SE,collapse=","),
                        sep=" "),""),
    				"2>",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"out",sep=".")),
    				"| samtools view -bS -F 0x100 - 2> /dev/null >", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),sep=" "));
    	})
    },mc.cores=floor(procs/opt$mprocs))
    if (!all(sapply(bowtie_out, "==", 0L))){
        write(paste("Something went wrong with bowtie mapping some jobs failed"),stderr())
        stop()
    }
}else if (opt$mappingAlgorithm == "bwa"){
    ## run bwa
    bwa_out <- mclapply(mapping, function(index){
        dir.create(file.path(opt$mappingFolder,index$sampleFolder),showWarnings=FALSE)
        res = res2 = res3 = 0        
        if(length(index$PE1)) try({
            res<-system(paste("bwa mem",
                         "-M", # Mark shorter split hits as secondary (for Picard compatibility).
                         "-t", opt$mprocs,
                         "-R", paste("'@RG",
                            paste("ID",index$sampleFolder,sep=":"),         
                            paste("SM",index$sampleFolder,sep=":"),				 
                            paste("PL","ILLUMINA",sep=":"),				 
                            paste("LB","whatever",sep=":"),				 
                            paste("PU","whatever",sep=":"),
                            paste("DS","Paired",sep=":"), "'",sep="\t"),			 
                         index$target_path,
                         paste(index$PE1,collapse=","),
                         paste(index$PE2,collapse=","),
                         "2>",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","out",sep=".")),
                         "| samtools view -bS -F 0x100 - 2> /dev/null >", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")),sep=" "));
            system(paste("samtools view  -H", 
                         file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")),
                         "| head -n -1 > ",
                         file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"header",sep=".")),
                         "2> /dev/null",sep=" "))
        })
        if(length(index$SE)) try({
            res2<-system(paste("bwa mem",
                         "-M", # Mark shorter split hits as secondary (for Picard compatibility).
                         "-t", opt$mprocs,
                         "-R", paste("'@RG",
                             paste("ID",index$sampleFolder,sep=":"),         
                             paste("SM",index$sampleFolder,sep=":"),    			 
                             paste("PL","ILLUMINA",sep=":"),				 
                             paste("LB","whatever",sep=":"),				 
                             paste("PU","whatever",sep=":"),
                             paste("DS","Paired",sep=":"),"'", sep="\t"),			 
                         index$target_path,
                         paste(index$SE,collapse=","),
                         "2>",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","out",sep=".")),
                         "| samtools view -bS -F 0x100 - 2> /dev/null >", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep=".")),sep=" "));
            system(paste("samtools view  -H", 
                         file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")),
                         "| head -n -1 > ",
                         file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"header",sep=".")),
                         "2> /dev/null",sep=" "))
        })
        ## reheader
        if (file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep="."))))
            file.remove(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")))
        if (length(index$PE1) & length(index$SE)){
            res3 <- system(paste("samtools cat -h", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"header",sep=".")), "-o",
                         file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
                         ifelse(file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep="."))),
                                file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")),""),
                         ifelse(file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep="."))),
                                file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep=".")),""),
                         "2> /dev/null",sep=" "))
        } else if (file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")))){
            res3 <- system(paste("samtools reheader", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"header",sep=".")),
                                 file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")),
                                 ">", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
                                 "2> /dev/null",sep=" "))
        } else if (file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep=".")))){
            res3 <- system(paste("samtools reheader", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"header",sep=".")),
                                 file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep=".")),
                           ">", file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
                           "2> /dev/null",sep=" "))
        } else {
            res3 = 1
        }
        if (file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep="."))))
            file.remove(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"PE","bam",sep=".")))     
        if (file.exists(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep="."))))
            file.remove(file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"SE","bam",sep=".")))     
        return(as.integer(res | res2 | res3))
    },mc.cores=floor(procs/opt$mprocs))
    if (!all(sapply(bwa_out, "==", 0L))){
        write(paste("Something went wrong with bwa mapping some jobs failed"),stderr())
        stop()
    }
}


## run samtools
samtools_out <- mclapply(mapping, function(index){
	dir.create(file.path(opt$mappingFolder,index$sampleFolder))
	try({
		res <- system(paste("samtools sort",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,sep=".")),"2> /dev/null",sep=" "));
		res <- res & system(paste("samtools index",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),"2> /dev/null",sep=" "));
		res <- res & system(paste("samtools idxstats",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),">",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"idxstats",sep=".")),"2> /dev/null",sep=" "))
		res <- res & system(paste("samtools flagstat",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),">",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"flagstat",sep=".")),"2> /dev/null",sep=" "))
		if (opt$sortByReadID) {
		    res <- res & system(paste("samtools sort -n",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"byreadid",sep=".")),"2> /dev/null",sep=" "));	
		}
        return(res)
	})
},mc.cores=procs)
if (!all(sapply(samtools_out, "==", FALSE))){
    write(paste("Something went wrong with samtools processing some jobs failed"),stderr())
    stop()
}

#####################################################
## write out index tables
targetTables <- sapply(targets,function(tgt){
	print(paste("Generating output for target:",tgt[1]))
	filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$mappingFolder,x,paste(x,tgt[1],"idxstats",sep="."))))
#	filesToRead <- unlist(sapply(file.path(opt$mappingFolder,unique(samples[,opt$samplesColumn])),dir,pattern=paste(tgt[1],"idxstats",sep="."),full.names=TRUE))
	info <- read.table(filesToRead[1])[,1:2]
	colnames(info) <- c("SequenceID","SequenceLength")
	data <- sapply(filesToRead,function(file){
		tb <- read.table(file)
		values <- rowSums(as.matrix(tb[,3:4]))
		values
	})
	colnames(data) <- basename(colnames(data))
	freq <- round(sweep(data,2,colSums(data),"/")*100,3)
	write.table(cbind(info,data),file.path(opt$mappingFolder,paste(tgt[1],"summary","reads","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	write.table(cbind(info,freq),file.path(opt$mappingFolder,paste(tgt[1],"summary","proportions","txt",sep=".")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	pmapped <- 100-tail(freq,1)
#	as.vector(pmapped)
	names(pmapped) <- colnames(freq)
	round(pmapped,3)
})

targetTables <- data.frame(targetTables)
colnames(targetTables) <- sapply(targets,"[[",1L)

### simple assign by most on target
targetTables <- data.frame(ID=rownames(targetTables),targetTables,assign=colnames(targetTables)[apply(targetTables,1,which.max)])
write.table(targetTables,file.path(opt$mappingFolder,"SummarySample2Targets.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

filesToRead <- unlist(sapply(unique(samples[,opt$samplesColumn]),function(x) file.path(opt$mappingFolder,x,paste(x,index$target_name,"flagstat",sep="."))))
data.flagstat <- sapply(filesToRead,function(file){
    values <- readLines(file)
    values <- sapply(strsplit(values,split=" + 0",fixed=T),"[[",1L)
    as.numeric(values)
})
rownames(data.flagstat) <- c("totalNumberOfReads","duplicates","numMappedReads","ReadsPaired","read1","read2","ProperlyPaired","itselfandmate","Singletons","mappedAcrossContigs","mapChrQ5")
data.flagstat = rbind(data.flagstat[c("totalNumberOfReads","numMappedReads"),],"numMappedReadsPercent"=data.flagstat["numMappedReads",]/data.flagstat["totalNumberOfReads",],data.flagstat[c("ReadsPaired","ProperlyPaired"),],"ProperlyPairedPercent"=data.flagstat["ProperlyPaired",]/data.flagstat["ReadsPaired",],data.flagstat[c("Singletons","mappedAcrossContigs"),])
write.table(t(data.flagstat),file.path(opt$mappingFolder,"MappingFlagstats.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)





#####################################################
## extract unmapped reads
if (opt$extract_unmapped){## Extract Unmapped Reads
	dir.create(file.path(opt$screenFolder))
	extract_out <- mclapply(mapping, function(index){
		try({
				dir.create(file.path(opt$screenFolder,index$sampleFolder));
				system(paste("samtools view",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_unmapped_reads.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"unmapped",sep=".")),sep=" "),intern=TRUE);
		})
	},mc.cores=procs/2)
	extract_out <- strsplit(sapply(extract_out,tail,n=1),split=": |,")
	extract_table <- data.frame(ID=names(mapping),Records=sapply(extract_out,"[[",2L),PE_pairs=sapply(extract_out,"[[",4L),SE_reads=sapply(extract_out,"[[",6L))
	write.table(extract_table,file.path(opt$screenFolder,"SummaryUnmapped.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

# #####################################################
# ## extract unmapped reads
if (opt$extract_mapped){## Extract Unmapped Reads
	 dir.create(file.path(opt$screenFolder))
	 extract_out <- mclapply(mapping, function(index){
		 try({
			 dir.create(file.path(opt$screenFolder,index$sampleFolder));
			 system(paste("samtools view",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")), "| extract_mapped_reads.py",ifelse(opt$gzip_extracted,"","-u"),"-v -o",file.path(opt$screenFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"mapped",sep=".")),sep=" "),intern=TRUE);
		 })
	 },mc.cores=procs/2)
	 extract_out <- strsplit(sapply(extract_out,tail,n=1),split=": |,")
	 extract_table <- data.frame(ID=names(mapping),Records=sapply(extract_out,"[[",2L),PE_pairs=sapply(extract_out,"[[",4L),SE_reads=sapply(extract_out,"[[",6L))
	 write.table(extract_table,file.path(opt$screenFolder,"SummaryMapped.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

######################################################
## generate VCF files

#java -Xmx6g -jar /mnt/home/msettles/opt/src/GenomeAnalysisTK-2.7-4/GenomeAnalysisTK.jar	-T HaplotypeCaller -R ../../03-targets/A.inornata/A.inornata_combined.fasta -I CP1.A.inornata.bam -o CP1.A.inornata.vcf
# 
# gtk <- "java -Xmx6g -jar /mnt/home/msettles/opt/src/GenomeAnalysisTK-2.7-4/GenomeAnalysisTK.jar"
# if (opt$generate_vcf){
#	 vcf_out <- mclapply(mapping, function(index){
#		 try({
#			 system(paste(gtk,
#									"-T HaplotypeCaller",
#									"-R",paste(index$target_path,"fasta",sep="."),
#									"-I",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"bam",sep=".")),
#									"-o",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"vcf",sep=".")),
#									">",file.path(opt$mappingFolder,index$sampleFolder,paste(index$sampleFolder,index$target_name,"gtk","out",sep=".")),sep=" "));		
#		 })
#	 },mc.cores=procs)
# }
