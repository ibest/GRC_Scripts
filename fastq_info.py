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
            c = 0
            seq1 = iterator1.next()
            count += 1
            lcount += 1
            bases += len(seq1)
            lbases += len(seq1)
            lqual += sum(read1.letter_annotations['phred_quality'])/lbases
            
    except StopIteration:
        pass
    finally:
        print "Finished processing file" +  infile
        outfile1.write("file\t" + infile1)
        outfile1.write("nreads\t" + lcount)
        outfile1.write("nbases\t" + lbases)
        outfile1.write("avgBases\t" + lbases/lcount)
        outfile1.wrtie("avgQual\t" + lqual/lcount)

#main part of the program

count = 0
bases = 0

outfile1 = open(os.path.realpath(os.path.join(os.getcwd(), sample_dir, "read_data.txt")))    
files = listdir_nohidden('./' + sample_dir)
    
for f in files:
    if ["fastq","fq"] in f
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        main(infile1, outfile1)

outfile1.write("directory\t" + sample_dir)
outfile1.write("treads\t" + count)
outfile1.write("tbases\t" + bases)


outfile1.close()

