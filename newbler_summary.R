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
    make_option(c("-n", "--newblerFolder"), type="character", default="03-NewblerAssemblies",
                help="Directory where to store the newbler results [default %default]",
                dest="newblerFolder")
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- list(samplesFile="samples.txt",readsFolder="02-Cleaned",newblerFolder="05-Assemblies")
#opt <- list(samplesFile="samples.txt", samplesColumn="SAMPLE_ID", readFolder="04-Screened-extra_human",newblerFolder="05-NewblerAssemblies",  procs=0, nprocs=10)

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
"loadSamplesFile" <- function(file,column){
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
    write(paste("samples sheet contains", nrow(targets), "samples to process",sep=" "),stdout())  
    return(targets)  
}


samples <- loadSamplesFile(opt$samplesFile,opt$samplesColumn)

"parse_newblerFiles" <- function(file){
    lines <- readLines(file)
    # remove extra whitespace
    lines <- gsub("\\s+$", "", lines)
    lines <- lines[!(lines == "")]
    # remove comments
    scomment <- grep("/*",lines,fixed=T)
    ecomment <- grep("*/",lines,fixed=T)
    if (length(scomment) != length(ecomment)){ 
        write(paste("Comment tags aren't matches, file:",file,"\n"), stderr())
        stop()
    }
    lines <- lines[-unlist(apply(cbind(scomment,ecomment),1,function(x) seq.int(x[1],x[2])))]
    
    procLines <- function(iloop){
        res <- {}
        slist <- which(iloop== "{")
        elist <- which(iloop == "}")
        if (length(slist) > 0 & length(slist) == length(elist)){
            nlist <- split(iloop[unlist(apply(cbind(slist,elist),1,function(x) seq.int(x[1]+1,x[2]-1)))],
                           rep(seq.int(1,length(slist)),times=apply(cbind(slist,elist),1,function(x) (x[2]-x[1]-1))))
            names(nlist) = iloop[slist-1]
            res <- lapply(nlist,function(x) procLines(gsub("^\t","",x)))
            iloop = iloop[-unlist(apply(cbind(slist,elist),1,function(x) seq.int(x[1]-1,x[2])))]
        }
        if (length(iloop) > 0){
            iloop <- gsub("^\\s+", "", iloop)
            iloop <- gsub(";|\"| MB", "", iloop)
            ll <- strsplit(iloop,split=" += +")
            nlist <- lapply(ll,"[[",2L)
            names(nlist) <- sapply(ll,"[[",1L)
            res <- c(nlist,res)
        }
        return(res)
    }
    
    # outer list  
    return(procLines(lines))
}

newblertb <- sapply(samples[,opt$samplesColumn], function(newb){
    require("Hmisc")
    cfile <- t(data.frame(strsplit(grep("contig",readLines(file.path(opt$newblerFolder,newb,"454ContigGraph.txt")),value=TRUE),"\t"),stringsAsFactors=F))
    cfile <- data.frame("Contig"=cfile[,2],"Length"=as.numeric(cfile[,3]),"Cov"=as.numeric(cfile[,4]))
    cov <- c(wtd.mean(cfile$Cov,cfile$Length),sqrt(wtd.var(cfile$Cov,cfile$Length)))    
    pfile <- parse_newblerFiles(file.path(opt$newblerFolder,newb,"454NewblerMetrics.txt"))
    # run data
    areads <- as.numeric(c(pfile$runMetrics$totalNumberOfReads, unlist(strsplit(gsub("%","",pfile$consensusResults$readStatus$numAlignedReads),split=" *, *"))))
    abases <- as.numeric(c(pfile$runMetrics$totalNumberOfBases, unlist(strsplit(gsub("%","",pfile$consensusResults$readStatus$numAlignedBases),split=" *, *"))))
    rstatus <- c(pfile$consensusResults$readStatus$numberAssembled,
                            pfile$consensusResults$readStatus$numberPartial,
                            pfile$consensusResults$readStatus$numberSingleton,
                            pfile$consensusResults$readStatus$numberRepeat,
                            pfile$consensusResults$readStatus$numberOutlier,
                            pfile$consensusResults$readStatus$numberTooShort)
    rstatus <- as.numeric(sapply(strsplit(rstatus,split=" |, "),"[[",1L))
    passembled <- (sum(rstatus[0:1])/areads[1])*100
    largecontigs <- as.numeric(c(pfile$consensusResults$largeContigMetrics$numberOfContigs,
                                 pfile$consensusResults$largeContigMetrics$numberOfBases,
                                 pfile$consensusResults$largeContigMetrics$avgContigSize,
                                 pfile$consensusResults$largeContigMetrics$N50ContigSize,
                                 pfile$consensusResults$largeContigMetrics$largestContigSize,
                                 unlist(strsplit(gsub("%","",pfile$consensusResults$largeContigMetrics$Q40PlusBases),split=" *, *"))))
    allcontigs <- as.numeric(c(pfile$consensusResults$allContigMetrics$numberOfContigs,pfile$consensusResults$allContigMetrics$numberOfBases))
    ndata <- c(areads[1:3],abases[1:3],rstatus,passembled,largecontigs,allcontigs,cov)
    names(ndata) <- c("totalNumberOfReads","numAlignedReads","numAlignedReadsPercent",
                      "totalNumberOfBases","numAlignedBases","numAlignedReadsBases",
                      "numberAssembled","numberPartial","numberSingleton","numberRepeat","numberOutlier","numberTooShart","assembledPercent",
                      "numLargeContigsAssembled","numLargeBasesAssembled","avgLargeContigSize","N50LargeContigSize","largestContigSize","numQ40PlusBases","Q40PlusBasesPercent",
                      "numAllContigsAssembled","numAllBasesAssembled","meanWeightedCov","sdWeightedCov")
    return(round(ndata,3))
})

newblertb <- t(newblertb)
write.table(newblertb,file.path(opt$newblerFolder,"SummaryNewblerAssemblies.txt"),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
