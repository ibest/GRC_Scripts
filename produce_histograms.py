#!/usr/bin/env python

"""
Very simply script to create a histogram from reads/or mapping

# Copyright 2016, Matt Settles
# Modified Jan 1, 2016
"""

import gzip
from Bio import SeqIO
import sys
import numpy

if len(sys.argv) != 4:
    print "Usage: produce_histograms.py file1.fastq.gz file2.fastq.gz file3_SE.fastq.gz"
    sys.exit()

fq1 = sys.argv[1]
fq2 = sys.argv[2]
fqse = sys.argv[3]

# fq1 = "BS_Eggs_R1.fastq"
# fq2 = "BS_Eggs_R2.fastq"
# fqse = "BS_Eggs_SE.fastq"

if fq1.split('.')[-1] == 'gz':
    h1 = SeqIO.parse(gzip.open(fq1, 'rb'), 'fastq')
else:
    h1 = SeqIO.parse(open(fq1, 'r'), 'fastq')

if fq2.split('.')[-1] == 'gz':
    h2 = SeqIO.parse(gzip.open(fq2, 'rb'), 'fastq')
else:
    h2 = SeqIO.parse(open(fq2, 'r'), 'fastq')

if fqse.split('.')[-1] == 'gz':
    s1 = SeqIO.parse(gzip.open(fqse, 'rb'), 'fastq')
else:
    s1 = SeqIO.parse(open(fqse, 'r'), 'fastq')


pairs_i = 0
pairs_v = []
try:
    print "Analyzing Pairs"
    while True:
        r1 = h1.next()
        r2 = h2.next()
        pairs_v.extend([len(r1.seq) + len(r2.seq)])
        pairs_i += 1
except StopIteration:
    pass
print "Finished processing %s pairs of reads." % pairs_i 

singles_i = 0
singles_v = []
try:
    print "Analyzing Singles"
    while True:
        se1 = s1.next()
        singles_v.extend([len(se1.seq)])
        singles_i += 1
except StopIteration:
    pass
print "Finished processing %s single reads." % singles_i 

h1.close()
h2.close()
s1.close()

print "Producing histogram files"
puni = numpy.unique(pairs_v)
phist, bin_edges = numpy.histogram(pairs_v, bins=numpy.append(puni,puni[-1]+1) - 0.1)

suni = numpy.unique(singles_v)
shist, bin_edges = numpy.histogram(singles_v, bins=numpy.append(suni,suni[-1]+1) - 0.1)

def write_hist_file(hist_file, hist, unique):
    fp = open(hist_file, 'w')
    for i in range(min(unique), max(unique)+1):
        try:
            value = hist[numpy.where(unique==i)][0]
        except (ValueError, IndexError):
            value = ""
        fp.write( "%s\t%s\n" % (i, value))
    fp.close()


def write_histogram_file(histogram_file, hist, unique):
    max_num_asterisks = float(72)
    scale = float(max_num_asterisks) / float(max(hist))
    fp = open(histogram_file, 'w')
    fp.write("# each asterisk represents approximately %s reads\n" % int(max(hist)/max_num_asterisks))
    for i in range(min(unique), max(unique)+1):
        try:
            num_asterisks = int(scale*hist[numpy.where(unique==i)][0])
        except (ValueError, IndexError):
            num_asterisks = 0
        fp.write( "%s\t%s\n" % (i, '*' * num_asterisks))
    fp.close()

write_hist_file(os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.hist"), phist, puni)

write_histogram_file(os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.histogram"), phist, puni)

write_hist_file(os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.hist"), shist, suni)

write_histogram_file(os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.histogram"), shist, suni)

