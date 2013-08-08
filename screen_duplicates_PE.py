#!/usr/bin/env python
"""
# Copyright 2013, Sam Hunter
"""


from Bio import SeqIO
from optparse import OptionParser
import sys
from collections import Counter
import time

## Parse options and setup ##
usage = "usage %prog read1.fastq read2.fastq output_base"
parser = OptionParser(usage=usage)
#parser.add_option('-a', '--adapterlength', help="length of adapter, controls max dovetail overlap",
#                  action="store", type="int", dest="adapterlength", default=30)

(options, args) = parser.parse_args()

#adapterlength = options.adapterlength
#minoverlap = 10  # minimum value at which at 90% overlap will be accepted (10bp overlap, 9/10 match)

if len(args) != 3:
    parser.print_help()
    sys.exit()

infile1 = args[0]
infile2 = args[1]

#Open inputs:
if infile1.split(".")[-1] == "gz":
    import gzip
    iter1 = SeqIO.parse(gzip.open(infile1, 'rb'), 'fastq')
    iter2 = SeqIO.parse(gzip.open(infile2, 'rb'), 'fastq')
elif infile1.split(".")[-1] == "fastq":
    iter1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
    iter2 = SeqIO.parse(open(infile2, 'r'), 'fastq')
else:
    iter1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
    iter2 = SeqIO.parse(open(infile2, 'r'), 'fastq')

pe1_outf = open(args[2] + "_nodup_PE1.fastq", 'w')
pe2_outf = open(args[2] + "_nodup_PE2.fastq", 'w')

#pe2_outf = gzip.open(args[2] + "_nodup_PE2.fastq.gz", 'wb')


def main():
    #main part of the program
    count = Counter()
    i = 0
    duplicates = 0
    rev = 0
    start = time.time()
    try:
        while 1:
            c = 0
            seq1 = iter1.next()
            seq2 = iter2.next()
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
                SeqIO.write(seq1, pe1_outf, "fastq")
                SeqIO.write(seq2, pe2_outf, "fastq")
            else:
                duplicates += 1
            i += 1
            if i % 10000 == 0:
                print "Pairs:", i, "Duplicates:", duplicates, "| fw:", duplicates-rev, "| rev", rev, "| percent:", 100.0*duplicates/i, "| reads/second:", i/(time.time() - start)

    except StopIteration:
        pass
    finally:
        print "Finished processing"
        print "Final:", i, "Duplicates:", duplicates, "| fw:", duplicates-rev, "| rev", rev, "| percent:", 100.0*duplicates/i, "| reads/second:", i/(time.time() - start)

main()
