#!/usr/bin/env python
"""
# Copyright 2013, Sam Hunter, Brice Sarver, Matt Settles
# Modified Aug 9, 2013
"""


from Bio import SeqIO
from optparse import OptionParser
import sys, os, os.path, time, gzip
from collections import Counter

## Parse options and setup ##
usage = "usage %prog -d [path to directory of raw reads] -o [path to output directory]"
#usage = "usage %prog read1.fastq read2.fastq output_base"
parser = OptionParser(usage=usage)

parser.add_option('-d', '--directory', help="Directory containing read files to de-duplicate",
    action="store", type="str", dest="sample_dir")

parser.add_option('-o', '--output', help="Directory to output de-duplicated reads",
    action="store", type="str", dest="output_dir")


(options, args) = parser.parse_args()

if len(args) != 0:
    parser.print_help()
    sys.exit()

sample_dir = options.sample_dir
output_dir = options.output_dir

#kindly provided by http://stackoverflow.com/questions/7099290/how-to-ignore-hidden-files-using-os-listdir-python
#glob.glob will list hidden files
#this replaces that functionality when hidden files exist, like in my reads from Berkeley
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f


def main(infile1, infile2, outfile1, outfile2):
#Set up the global variables
    global count
    global i
    global duplicates
    global rev
    global start
#Open inputs:
    if infile1.split(".")[-1] == "gz":
        import gzip
        iterator1 = SeqIO.parse(gzip.open(infile1, 'rb'), 'fastq')
        iterator2 = SeqIO.parse(gzip.open(infile2, 'rb'), 'fastq')
    elif infile1.split(".")[-1] == "fastq":
        iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
        iterator2 = SeqIO.parse(open(infile2, 'r'), 'fastq')
    else:
        iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
        iterator2 = SeqIO.parse(open(infile2, 'r'), 'fastq')

    try:
        while 1:
            c = 0
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            #comb = seq1.seq.tostring() + seq2.seq.tostring()
            comb = seq1[10:35] + seq2[10:35]
            rcomb = comb.reverse_complement()
            comb = comb.seq.tostring()
            rcomb = rcomb.seq.tostring()
            if rcomb in count:
                count[rcomb] += 1
                c = count[rcomb]
                rev += 1
            else:
                count[comb] += 1
                c = count[comb]
            if c == 1:
                seq1.id = seq1.id + ":0:0:0:0:0:0 1:Y:0:"
                SeqIO.write(seq1, outfile1, "fastq")
                seq2.id = seq2.id + ":0:0:0:0:0:0 2:Y:0:"
                SeqIO.write(seq2, outfile2, "fastq")
            else:
                duplicates += 1
            i += 1
            if i % 100000 == 0:
              print "Pairs:","| reads:", i, "| duplicates:", duplicates, "| fwd:", duplicates-rev, "| rev:", rev, "| percent:", round(100.0*duplicates/i,2), "| reads/sec:", round(i/(time.time() - start),0)
            
    except StopIteration:
        pass
    finally:
#Print final output
        print "Finished processing one pair of files."

    #main part of the program
count = Counter()
i = 0
duplicates = 0
rev = 0
start = time.time()

outfile1 = gzip.open(output_dir + "_nodup_PE1.fastq.gz", 'wb')
outfile2 = gzip.open(output_dir + "_nodup_PE2.fastq.gz", 'wb')
    
files = listdir_nohidden('./' + sample_dir)
    
for f in files:
    if "_R1" in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "_R2".join(f.split("_R1"))))
        main(infile1, infile2, outfile1, outfile2)
    elif "_1.fastq" in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "_2.fastq".join(f.split("_1.fastq"))))
        main(infile1, infile2, outfile1, outfile2)

outfile1.close()
outfile2.close()

print "Final:","| reads:", i, "| duplicates:", duplicates, "| fwd:", duplicates-rev, "| rev:", rev, "| percent:", round(100.0*duplicates/i,2), "| reads/sec:", round(i/(time.time() - start),0)

