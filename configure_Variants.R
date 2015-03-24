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
    make_option(c("-r", "--variantFolder"), type="character", default="04-Variants",
                help="Directory where the variant data will be stored [default %default]",
                dest="variantFolder"),
    make_option(c("-M", "--mappingAlgorithm"), type="character", default="bwa",
                help="Mapping algorithm to use, supported types are 'bowtie' and 'bwa' [default %default]",
                dest="mappingAlgorithm"),
    make_option(c("-b", "--mappingFolder"), type="character", default=NA,
                help="Directory where to store the mapping results [default '03-[mappingAlgorithm]']",
                dest="mappingFolder"),
    make_option(c("-t", "--mappingTargets"), type="character", default="mapping_targets.txt",
                help="Path to a fasta file, or tab delimeted file with [target name]\t[target fasta]\t[target gtf, optional] to run mapping against [default %default]",
                dest="mappingTarget"),
#     make_option(c("-p", "--processors"), type="integer", default=0,
#                 help="number of processors to use [defaults to number available]",
#                 dest="procs"),
    make_option(c("-C", "--combined"), action="store_true", default=FALSE,
                help="combined the bam files for analysis, produce a single vcf file [default %default]",
                dest="combined"),
    make_option(c("--PICARDpath"), type="character", default="/mnt/home/grcuser/module_grc/src/PICARD/picard-tools-1.119/CreateSequenceDictionary.jar",
                help="Path to the PICARD jar files [default %default]",
                dest="PICARD"),
        make_option(c("--GATKpath"), type="character", default="/mnt/home/grcuser/module_grc/src/GATK/GenomeAnalysisTK.jar",
                help="Path to the GATK jar file [default %default]",
                dest="GATK")
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
"loadSamplesFile" <- function(file,column){
    ## debug
#    file = opt$samplesFile; reads_folder = opt$readFolder; column = opt$samplesColumn
    ##
    if ( !file.exists(file) ) {
        write(paste("Sample file",file,"does not exist\n"), stderr())
        stop()
    }	
    ### column SEQUENCE_ID should be the folder name inside of Raw_Folder
    ### column SAMPLE_ID should be the sample name
    ### rows can be commented out with #
    targets <- read.table(file,sep="",header=TRUE,as.is=TRUE)
    if( !all(column %in% colnames(targets)) ){
        write(paste("Expecting the two columns", paste(column,collapse=","),"in samples file (tab-delimited)\n"), stderr())
        stop()
    }
     return(targets)	
}

######################################################################
## prepareTargets
##	Prepare the mapping targets 
## 
## Parameters
##	targets: filename of targets builds, fasta file or text file with multiple targets
picard = paste("java -jar",opt$PICARD)
"prepareTargets" <- function(targets){
    ### single target, indexes exist
    if (file.exists(targets) & file_ext(targets) %in% c("fasta","fa","fna")){
        targets_list <- list(c(basename(targets),targets))
    } else if (file.exists(targets)){
        ### multiple targets
        targets_list <- lapply(readLines(targets),function(x) strsplit(x,split="\t")[[1]])
        #	 Assume first column is name, second is the fasta file, remaining columns are ignored
    } else {
        write(paste("Something wrong with targets file (or table)"),stderr())
        stop("Quiting")
    }
    write(paste("Found", length(targets_list), "targets",sep=" "),stdout())	
    sapply(targets_list,function(targets){
        if(!file.exists(paste(targets[2],"fai",sep=".")))
            system(paste("samtools faidx",targets[2]))
        if(!file.exists(ifelse(file_ext(targets[2]) %in% c("fasta","fa"),paste(sub(".fasta$|.fa$","",targets[2]),"dict",sep="."),paste(targets[2],"dict",sep="."))))
            system(paste(picard,"R=",targets[2],"O=",ifelse(file_ext(targets[2]) %in% c("fasta","fa"),paste(sub(".fasta$|fa$","",targets[2]),"dict",sep="."),paste(targets[2],"dict",sep="."))))        
    })    
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
"mappingList" <- function(samples, mapping_folder, column, targets, algorithm){
    mapping_list <- list()
    for (i in seq.int(to=nrow(samples))){
        map = list()
        map$sampleName=samples[i,column]
        for(j in targets){
            mapping_list[[paste(map$sampleName,j[1],sep="_")]] <- map
            mapping_list[[paste(map$sampleName,j[1],sep="_")]]$target_name <- j[1]
            if(opt$mappingAlgorithm == "bowtie")
                mapping_list[[paste(map$sampleName,j[1],sep="_")]]$target_name <- sub(".fasta$|.fa$|.fna$","",j[1])
            mapping_list[[paste(map$sampleName,j[1],sep="_")]]$target_path <- j[2]
            mapping_list[[paste(map$sampleName,j[1],sep="_")]]$bam_file <- file.path(mapping_folder,map$sampleName,paste(map$sampleName,j[1],"bam",sep="."))
        }
    }
    write(paste("Setting up",length(mapping_list),"jobs",sep=" "),stdout())
    return(mapping_list)
}


samples <- loadSamplesFile(opt$samplesFile,opt$samplesColumn)
targets <- prepareTargets(opt$mappingTarget)
#procs <- prepareCore(opt$procs)
mapping <- mappingList(samples,opt$mappingFolder,opt$samplesColumn,targets,opt$mappingAlgorithm)

## create output folder
out <- sapply(file.path(opt$variantFolder,sapply(mapping,"[[","sampleName")),dir.create,showWarnings=FALSE,recursive=TRUE)

######################################################
## generate VCF files
gtk <- paste("java -Xmx6g -jar",opt$GATK)


indelrealign <- mclapply(mapping,function(index){
    try({
        call <- paste(gtk,
                      "-T RealignerTargetCreator",
                      "-R", index$target_path,
                      "-I", index$bam_file,
                      "-o", file.path(opt$variantFolder,index$sampleName,paste(basename(index$bam_file),"intervalListFromRTC.intervals",sep=".")))
        system(call)
    })
},mc.cores=8)



indelrealign <- mclapply(mapping,function(index){
    try({
        call <- paste(gtk,
                      "-T IndelRealigner",
                      "-R", index$target_path,
                      "-I", index$bam_file,
                      "-targetIntervals", file.path(opt$variantFolder,index$sampleName,paste(basename(index$bam_file),"intervalListFromRTC.intervals",sep=".")),
                      "-o", file.path(opt$variantFolder,index$sampleName,paste(basename(index$bam_file),"realigned.bam",sep=".")))
        system(call)
    })
},mc.cores=8)


if (opt$combined){
        target_names <- sapply(mapping,"[[","target_name")
        mclapply(split(mapping,target_names),function(target_map){
            call <- paste(gtk,
                          "-T HaplotypeCaller",
                          "-R",target_map[[1]]$target_path,
                          paste(paste("-I",file.path(opt$variantFolder,sapply(target_map,"[[","sampleName"),paste(basename(sapply(target_map,"[[","bam_file")),"realigned.bam",sep="."))),collapse=" "),
                          "-o",file.path(opt$variantFolder,paste("combined",target_map[[1]]$target_name,"vcf",sep=".")),
                          ">",file.path(opt$variantFolder,paste("combined",target_map[[1]]$target_name,"gtk","out",sep=".")),sep=" ")
            system(call)            
        },mc.cores=8)
    
} else {
    vcf_out <- mclapply(mapping, function(index){
        try({
            call <- paste(gtk,
                          "-T HaplotypeCaller",
                          "-R",index$target_path,
                          "-I",file.path(opt$variantFolder,index$sampleName,paste(basename(index$bam_file),"realigned.bam",sep=".")),
                          "-o",file.path(opt$variantFolder,index$sampleName,paste(index$sampleName,index$target_name,"vcf",sep=".")),
                          ">",file.path(opt$variantFolder,index$sampleName,paste(index$sampleName,index$target_name,"gtk","out",sep=".")),sep=" ")
            system(call)
        })
    },mc.cores=8)   
}
