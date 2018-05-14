### preproc_experiment.R
# updated script for preprocessing HTS data


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
  make_option(c("-d", "--directory"), type="character", default="00-RawData",
              help="Directory where the raw sequence data is stored [default %default]",
              dest="Raw_Folder"),
  make_option(c("-q", "--quality"), type="integer", default=24,
              help="Quality score to use during lucy trimming [default %default]",
              dest="qual"),
  make_option(c("-m", "--miniumumLength"), type="integer", default=150,
              help="Discard reads less then minimum length [default %default]",
              dest="minL"),
  make_option(c("-o", "--overlap"), type="integer", default=275,
              help="Overlap parameter for flash [default %default]",
              dest="overlap"),
  make_option(c("-p", "--processors"), type="integer", default=1,
              help="number of processors to use [default %default]",
              dest="procs"),
  make_option(c("-s", "--skip-duduplicates"), action="store_true", default=FALSE,
              help="do not perform the deduplication step [default %default] NOT FUNCTIONAL YET",
              dest="skip_dedup"),
  make_option(c("-c", "--contaminants-folder"), type="character", default=NULL,
              help="folder name with contaminant sequences in fasta format [default %default]",
              dest="contaminants"),
  make_option(c("-v", "--vector-folder"), type="character", default=NULL,
              help="folder name with vector sequences in fasta format [default %default]",
              dest="vector"),
  make_option(c("--i64"), action="store_true",default=FALSE,
              help="input read Q scores are offset by 64 [default %default]",
              dest="i64")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

