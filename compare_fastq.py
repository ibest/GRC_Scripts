#!/usr/bin/env python

"""
Very simply script which compares two fastq files

# Copyright 2014, Sam Hunter
# Modified Aug 27, 2014
"""


import gzip
from Bio import SeqIO
import sys
import numpy

if len(sys.argv) != 3:
    print "Usage: compare_fastq.py file1.fastq.gz file2.fastq.gz"
    sys.exit()

fq1 = sys.argv[1]
fq2 = sys.argv[2]
fqse = sys.argv[3]

fq1 = "BS_Eggs_R1.fastq"
fq2 = "BS_Eggs_R2.fastq"
fqse = "BS_Eggs_SE.fastq"
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

pairs_hist = os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.hist")

pairs_histogram = os.path.join(os.path.dirname(os.path.abspath(fq1)),"pairs.histogram")

singles_hist = os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.hist")

singles_histogram = os.path.join(os.path.dirname(os.path.abspath(fqse)),"singles.histogram")

def write_hist_file(hist_file, hist, unique):
    fp = open(hist_file, 'w')
    for i in range(min(unique), max(unique)+1):
        fp.write( "%s\t%s\n" % (i, hist[numpy.where(unique==i)][0]))
    fp.close()

static void
write_histogram_file(const char *histogram_file, const struct histogram *hist,
                     long first_nonzero_idx, long last_nonzero_idx,
                     uint64_t max_freq)
{
        const double max_num_asterisks = 72;
        double scale = max_num_asterisks / (double)max_freq;

        FILE *fp = xfopen(histogram_file, "w");

        for (long i = first_nonzero_idx; i <= last_nonzero_idx; i++) {
                if (fprintf(fp, "%ld\t", i) < 0)
                        goto write_error;
                size_t num_asterisks = (size_t)(scale * (double)hist_count_at(hist, i));
                while (num_asterisks--)
                        if (fputc('*', fp) == EOF)
                                goto write_error;
                if (fputc('\n', fp) == EOF)
                        goto write_error;
        }
        xfclose(fp, histogram_file);
        return;
write_error:
        fatal_error_with_errno("Error writing to \"%s\"", histogram_file);
}

