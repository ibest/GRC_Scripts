#!/usr/bin/env python

"""
Very simply script which compares two fastq files

# Copyright 2014, Sam Hunter
# Modified Aug 27, 2014
"""


import gzip
from Bio import SeqIO
import sys

if len(sys) != 3:
    print "Usage: compare_fastq.py file1.fastq.gz file2.fastq.gz"
    sys.exit()

fq1 = sys.argv[1]
fq2 = sys.argv[2]

#fq1 = "./pre/test_nodup_PE1.fastq"
#fq2 = "./pre/test_nodup_PE2.fastq"

if fq1.split('.')[-1] == 'gz':
    h1 = SeqIO.parse(gzip.open(fq1, 'rb'), 'fastq')
else:
    h1 = SeqIO.parse(open(fq1, 'r'), 'fastq')

if fq2.split('.')[-1] == 'gz':
    h2 = SeqIO.parse(gzip.open(fq2, 'rb'), 'fastq')
else:
    h2 = SeqIO.parse(open(fq2, 'r'), 'fastq')


i = 0
same = 0
different = 0
try:
    while True:
        i += 1
        r1_1 = h1.next()
        r1_2 = h2.next()
        if str(r1_1.seq) != str(r1_2.seq):
            print "Mismatch record %s:" % i
            print '\t' + str(r1_1.seq)
            print '\t' + str(r1_2.seq)
            different += 1
        else:
            same += 1
except StopIteration:
    pass
finally:
    #Print final output
    print "Finished processing one pair of files."
    print "Total records: %s \n\tidentical records: %s  \n\tdifferent records: %s" % (i, same, different)
