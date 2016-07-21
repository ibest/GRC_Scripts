#!/usr/bin/env python
"""
# Copyright 2015 Matt Settles
# Created Dec 1, 2014
"""
from optparse import OptionParser
import sys
import time
from subprocess import Popen, PIPE, STDOUT
import string
import re
import copy


def sp_gzip_read(file, bufsize=-1):
    p = Popen('gzip --decompress --to-stdout'.split() + [file], stdout=PIPE, stderr=STDOUT, bufsize=bufsize)
    return p.stdout


def sp_gzip_write(file, bufsize=-1):
    filep = open(file, 'wb')
    p = Popen('gzip', stdin=PIPE, stdout=filep, shell=True, bufsize=bufsize)
    return p.stdin

rcs = string.maketrans('TAGCtagc', 'ATCGATCG')


def revcomp(seq):
    return seq.translate(rcs)[::-1]


def rev(seq):
    return seq[::-1]


class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle

    def __iter__(self):
        return self

    def next(self):
        lines = {'id': self.inf.readline().strip(),
                 'seq': self.inf.readline().strip(),
                 '+': self.inf.readline().strip(),
                 'qual': self.inf.readline().strip()}
        assert(len(lines['seq']) == len(lines['qual']))
        if lines['id'] == '' or lines['seq'] == '' or lines['+'] == '' or lines['qual'] == '':
            raise StopIteration
        else:
            return lines

    @staticmethod
    def parse(handle):
        return fastqIter(handle)

    def close(self):
        self.inf.close()


def writeFastq(handle, fq):
    handle.write(fq['id'] + '\n')
    handle.write(fq['seq'] + '\n')
    handle.write(fq['+'] + '\n')
    handle.write(fq['qual'] + '\n')


def main(read1, read2, single, outfile1, outfile2, min, length, verbose):
    # Set up the global variables
    global pair_count
    global pair_count_kept
    global pair_bases_count
    global pair_bases_count_kept
    global singles_count
    global singles_count_kept
    global singles_bases_count
    global singles_bases_count_kept
    global stime
    # Process Paired inputs:
    if read1.split(".")[-1] == "gz":
        iterator1 = fastqIter(sp_gzip_read(read1))
        # assume both are gz
        iterator2 = fastqIter(sp_gzip_read(read2))
    else:
        iterator1 = fastqIter(open(read1, 'r'))
        iterator2 = fastqIter(open(read2, 'r'))
    try:
        while 1:
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            pair_count += 1
            pair_bases_count += len(seq1['seq']) + len(seq1['seq'])
            if len(seq1['seq']) >= min and len(seq2['seq']) >= min:  # both reads must be > length
                pair_count_kept += 1
                pair_bases_count_kept += len(seq1['seq']) + len(seq1['seq'])
                writeFastq(outfile1, seq1)
                writeFastq(outfile2, seq2)
            if pair_count % 500000 == 0 and verbose:
                print "Pairs:", "| reads:", pair_count, "| kept:", pair_count_kept, "| percent:", round(100.0 * pair_count_kept / pair_count, 2), "| reads/sec:", round(pair_count / (time.time() - stime), 0)
    except StopIteration:
        if verbose:
            print "Pairs:", "| reads:", pair_count, "| kept:", pair_count_kept, "| percent:", round(100.0 * pair_count_kept / pair_count, 2), "| reads/sec:", round(pair_count / (time.time() - stime), 0)
        pass
    # Process Single Reads
    if singles.split(".")[-1] == "gz":
        iterator3 = fastqIter(sp_gzip_read(singles))
    else:
        iterator3 = fastqIter(open(singles, 'r'))
    try:
        while 1:
            seq3 = iterator3.next()
            singles_count += 1
            singles_bases_count += len(seq3['seq'])
            if len(seq3['seq']) <= length and len(seq3['seq']) >= min:
                singles_count_kept += 1
                singles_bases_count_kept += 2 * len(seq3['seq'])
                writeFastq(outfile1, seq3)
                seq4 = copy.copy(seq3)
                seq4['id'] = re.sub(' 1', ' 2', seq3['id'])
                seq4['seq'] = revcomp(seq3['seq'])
                seq4['qual'] = rev(seq3['qual'])
                writeFastq(outfile2, seq4)
            elif len(seq3['seq']) > length:
                singles_count_kept += 1
                singles_bases_count_kept += 2 * length
                seq3['seq'] = seq3['seq'][0:length]
                seq3['qual'] = seq3['qual'][0:length]
                writeFastq(outfile1, seq3)
                seq4 = copy.copy(seq3)
                seq4['id'] = re.sub(' 1', ' 2', seq3['id'])
                seq4['seq'] = revcomp(seq3['seq'])[0:length]
                seq4['qual'] = rev(seq3['qual'])[0:length]
                writeFastq(outfile2, seq4)
            if singles_count % 500000 == 0 and verbose:
                print "Singles:", "| reads:", singles_count, "| kept:", singles_count_kept, "| percent:", round(100.0 * singles_count_kept / singles_count, 2), "| reads/sec:", round(singles_count / (time.time() - stime), 0)
    except StopIteration:
        if verbose:
            print "Singles:", "| reads:", singles_count, "| kept:", singles_count_kept, "| percent:", round(100.0 * singles_count_kept / singles_count, 2), "| reads/sec:", round(singles_count / (time.time() - stime), 0)
        pass


#####################################
# Parse options and setup #
usage = "usage %prog -o [output file prefix (path + name)] -(bl) --quite -1 [read1] -2 [read2] -U [singles]"
usage += "PostFlash2Pairs.py will process read file produced by flash and generate paired end only reads"
parser = OptionParser(usage=usage, version="%prog 0.0.1")

parser.add_option('-o', '--output', help="Directory + prefix to output reads",
                  action="store", type="str", dest="output_dir", default="paired-reads")

parser.add_option('-b', '', help="minimum length to output",
                  action="store", type="int", dest="minimum", default=150)

parser.add_option('-l', '--length', help="original length of read",
                  action="store", type="int", dest="length", default=250)

parser.add_option('-1', '--Read1', help="Read1, check length and ouput",
                  action="store", type="str", dest="read1", default=None)

parser.add_option('-2', '--Read2', help="Read2, check length and output",
                  action="store", type="str", dest="read2", default=None)

parser.add_option('-U', '--Single', help="Single Reads, split into pairs and output",
                  action="store", type="str", dest="singles", default=None)

parser.add_option('-g', '--gzip', help="gzip the output",
                  action="store_true", dest="gzip", default=False)

parser.add_option('--quite', help="turn off verbose output",
                  action="store_false", dest="verbose", default=True)

(options, args) = parser.parse_args()

output_dir = options.output_dir
minimum = options.minimum
length = options.length

infile1 = options.read1
if infile1 is None:
    sys.stdout.write("Paired end file 1 is missing\n")
    sys.exit(1)
infile2 = options.read2
if infile2 is None:
    sys.stdout.write("Paired end file 2 is missing\n")
    sys.exit(1)
singles = options.singles
if singles is None:
    sys.stdout.write("Singles file is missing\n")
    sys.exit(1)

verbose = options.verbose


pair_count = 0
pair_count_kept = 0
pair_bases_count = 0
pair_bases_count_kept = 0

singles_count = 0
singles_count_kept = 0
singles_bases_count = 0
singles_bases_count_kept = 0

stime = time.time()

if options.gzip:
    outfile1 = sp_gzip_write(output_dir + "_PE1.fastq.gz")
    outfile2 = sp_gzip_write(output_dir + "_PE2.fastq.gz")
else:
    outfile1 = open(output_dir + "_PE1.fastq", "w")
    outfile2 = open(output_dir + "_PE2.fastq", "w")

main(infile1, infile2, singles, outfile1, outfile2, minimum, length, verbose)

outfile1.close()
outfile2.close()

print "Final:", "| total reads input:", singles_count + 2 * pair_count, "| total bases input:", singles_bases_count + pair_bases_count, "| total reads output:", 2 * singles_count_kept + 2 * pair_count_kept, "| total bases output:", singles_bases_count_kept + pair_bases_count_kept

sys.exit(0)
