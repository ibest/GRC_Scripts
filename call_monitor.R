#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list <- list(
  make_option(c("-c", "--call"), type="character", default=NULL,
              help="System Call to Profile",
              dest="call", metavar="CALL"),
  make_option(c("-o", "--output"), type="character", default="profile.txt",
              help="output file name to send profile data to, default is profile.txt",
              dest="output", metavar="FILENAME")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# test
# opt <- list(call="preproc_experiment -p 4",output="profile.txt")

my_pid <- Sys.getpid()
my_user <- Sys.info()["user"]

p <- mcparallel({system(paste(opt$call,"> log.txt"),ignore.stdout=TRUE,wait=TRUE)})

childPID <- function(pid) as.numeric(system(paste("pgrep -P",my_pid),intern=TRUE))
parallel:::mckill(as.numeric(system(paste("pgrep -P",my_pid),intern=TRUE)))


time_data = data.frame()

while(TRUE){
  if(user == ""){
    system('ps xo euser,rss,pid,args | grep ARC | grep python | grep -v "grep" > memuse.tsv')
  }else{
    system(paste('ps xo euser,rss,pid,args | grep ARC | grep python | grep', user,  '| grep -v "grep" > memuse.tsv'))
  }
  if(file.exists("memuse.tsv")){
    dat = read.table("memuse.tsv")
    colnames(dat) = c("User","Mem","PID","ARGS")
    dat = cbind(dat, time = rep(proc.time()[3], dim(dat)[1]))
    time_data = rbind(time_data, dat[,c("Mem","PID","time")])
    write.table(time_data, "time_data.tsv")
  }
  Sys.sleep(10)
}