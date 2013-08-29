#!/usr/bin/env python
"""
# Copyright 2013, Matt Settles
# Modified Aug 10, 2013
"""


from Bio import SeqIO
from optparse import OptionParser
import sys, os, os.path, time, gzip
from collections import Counter

## Parse options and setup ##
usage = "usage %prog -d [path to directory of raw reads]"
parser = OptionParser(usage=usage)

parser.add_option('-d', '--directory', help="Directory containing read files to de-duplicate",
    action="store", type="str", dest="sample_dir")

(options, args) = parser.parse_args()

if len(args) != 0:
    parser.print_help()
    sys.exit()

sample_dir = options.sample_dir
output_dir = options.sample_dir

#kindly provided by http://stackoverflow.com/questions/7099290/how-to-ignore-hidden-files-using-os-listdir-python
#glob.glob will list hidden files
#this replaces that functionality when hidden files exist, like in my reads from Berkeley
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f

def main_sff(infile, outfile1):
# Set up the global variables
    global count
    global bases
    lcount = 0
    lbases = 0
    lqual = 0

# Open inputs:
    iterator1 = SeqIO.parse(open(infile, 'r'), "sff")
    try:
        while 1:
            seq1 = iterator1.next()
            count += 1
            lcount += 1
            len = seq1.annotations["clip_qual_right"] - seq1.annotations["clip_qual_left"]
            bases += len
            lbases += len
            lqual += sum(seq1.letter_annotations['phred_quality'][seq1.annotations["clip_qual_left"]:seq1.annotations["clip_qual_right"]])/len

    except StopIteration:
        pass
    finally:
        print "Finished processing file" +  infile1
        outfile1.write("file\t" + infile1 + "\n")
        outfile1.write("nreads\t" + str(lcount) + "\n")
        outfile1.write("nbases\t" + str(lbases) + "\n")
        outfile1.write("avgBases\t" + str(round(lbases/lcount,0)) + "\n")
        outfile1.write("avgQual\t" + str(round(lqual/lcount,1)) + "\n")

def main(infile1, outfile1):
#Set up the global variables
    global count
    global bases
    lcount = 0
    lbases = 0
    lqual = 0
    
#Open inputs:
    if infile1.split(".")[-1] == "gz":
        import gzip
        iterator1 = SeqIO.parse(gzip.open(infile1, 'rb'), 'fastq')
    elif infile1.split(".")[-1] == "fastq":
        iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
    else:
        iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')

    try:
        while 1:
            seq1 = iterator1.next()
            count += 1
            lcount += 1
            bases += len(seq1)
            lbases += len(seq1)
            lqual += sum(seq1.letter_annotations['phred_quality'])/len(seq1)
            
    except StopIteration:
        pass
    finally:
        print "Finished processing file" +  infile1
        outfile1.write("file\t" + infile1 + "\n")
        outfile1.write("nreads\t" + str(lcount) + "\n")
        outfile1.write("nbases\t" + str(lbases) + "\n")
        outfile1.write("avgBases\t" + str(round(lbases/lcount,0)) + "\n")
        outfile1.write("avgQual\t" + str(round(lqual/lcount,1)) + "\n")

#main part of the program

count = 0
bases = 0

outfile1 = open(os.path.realpath(os.path.join(os.getcwd(), sample_dir, "read_data.txt")),"w+")
files = listdir_nohidden('./' + sample_dir)
    
for f in files:
    if ("fastq" in f) or ("fq" in f):
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        main(infile1, outfile1)
    if ("sff" in f):
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        main_sff(infile1, outfile1)

outfile1.write("directory\t" + sample_dir + "\n")
outfile1.write("treads\t" + str(count) + "\n")
outfile1.write("tbases\t" + str(bases) + "\n")


outfile1.close()

