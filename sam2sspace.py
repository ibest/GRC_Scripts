#!/usr/bin/env python
#
# MATT SETTLES
# SSPACE_sam-bam2tab.py

import sys
import os

from subprocess import Popen
from subprocess import PIPE
import signal
import shlex

import re
import warnings

'''
NOTES
1 QNAME String
2 FLAG Int
3 RNAME String
4 POS Int
5 MAPQ Int
6 CIGAR String
7 RNEXT String
8 PNEXT Int
9 TLEN Int
10 SEQ String
11 QUAL String

Note (the & operation produces a value equivalent to the binary with a 1 in that position):
0x1 template having multiple segments in sequencing
0x2 each segment properly aligned according to the aligner
0x4 segment unmapped
0x8 next segment in the template unmapped
0x10 SEQ being reverse complemented
0x20 SEQ of the next segment in the template being reversed
0x40 the first segment in the template
0x80 the last segment in the template
0x100 secondary alignment
0x200 not passing quality controls
0x400 PCR or optical duplicate

SSPACE compatable Tab format
  -contig of read 1
  -start position of read 1
  -end position of read 1
  -contig of read 2
  -start position of read 2
  -end position of read 2
'''


def sp_bam2sam(file):
    p = Popen(shlex.split('samtools view') + [file],
              stdout=PIPE,
              stderr=None,
              bufsize=-1,
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


def write_pair_tab(r1, r2):
    sep = '\t'
    # output.write(sep.join(map(str, r1[1:] + r2[1:])) + '\t' + r1[0] + r2[0] + '\n')
    output.write(sep.join(map(str, r1[1:] + r2[1:])) + '\n')


class cigarString:
    """ Class to parse and handle sam formatted cigar strings """
    pattern = re.compile('([MIDNSHPX=])')

    def __init__(self, cigar):
        values = self.pattern.split(cigar)[:-1]
        self.paired = (values[n:n + 2] for n in xrange(0, len(values), 2))  # pair values by twos

    def getAlignmentLength(self):
        g = 0
        for pair in self.paired:
            l = int(pair[0])
            t = pair[1]
            if t == 'M':
                g += l
            elif t == 'I':
                pass
            elif t == 'D':
                g += l
            elif t == 'N':
                pass
            elif t == 'S':
                pass
            elif t == 'H':
                pass
            elif t == 'P':
                pass
            elif t == 'X':
                pass
            elif t == '=':
                pass
            else:
                warnings.warn("encountered unhandled CIGAR character %s" % t)
                pass
        return g


if len(sys.argv) == 2:
    infile = sys.argv[1]
    if not os.path.exists(infile):
        print "Error, can't find input file %s" % infile
        sys.exit(1)

    if infile.split(".")[-1] == "bam":
        insam = sp_bam2sam(infile)
    elif infile.split(".")[-1] == "sam":
        insam = open(infile, 'r')
    else:
        print("Error, requires a sam/bam (ends in either .sam or .bam) file as input")
        sys.exit(1)
    filen, ext = os.path.splitext(infile)
    output = open(filen + ".tab", 'w')
else:
    # reading/writing from stdin/stdout
    insam = sys.stdin
    output = sys.stdout

i = 0
secondary_alignment_count = 0
single_end_count = 0
paired_end_count = 0
paired_end_alignment_count = 0
paired_end_broken = 0

PE1 = {}
PE2 = {}

for line in insam:
    # Comment/header lines start with @
    if i % 1000000 == 0 and i > 0:
        print "Records processed: %s, discarded secondary alignments: %s, single end reads: %s, paired end reads: %s, broken pairs: %s" % (i / 2, secondary_alignment_count, single_end_count, paired_end_count, paired_end_broken)

    if line[0] != "@" and len(line.strip().split()) >= 11:  # ignore header lines
        i += 1
        line2 = line.strip().split()
        flag = int(line2[1])

        if (flag & 0x100):  # secondary alignment
            secondary_alignment_count += 1
            continue
        # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
        # Handle SE:
        if not (flag & 0x1) and not (flag & 0x4):  # If a single end read present, ignore it
            single_end_count += 1
            continue
        # Handle PE:
        # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
        if (flag & 0x1) and ((flag & 0x4) or (flag & 0x8)):
            # one of the two reads in the pair is unmapped
            if (flag & 0x40):
                paired_end_count += 1
                paired_end_broken += 1
            if (flag & 0x4):
                paired_end_alignment_count += 1
            continue
        # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
        if (flag & 0x1) and not (flag & 0x4) and not (flag & 0x8):
            paired_end_alignment_count += 1
            ID = line2[0].split("#")[0]
            cigar = cigarString(line2[5])
            if (flag & 0x40):  # If read 1
                paired_end_count += 1
                if (flag & 0x10):  # if RevComp
                    r1 = ['R', line2[2], int(line2[3]) + cigar.getAlignmentLength() - 1, line2[3]]
                else:
                    r1 = ['F', line2[2], line2[3], int(line2[3]) + cigar.getAlignmentLength() - 1]
                if ID in PE2:
                    write_pair_tab(r1, PE2[ID])
                    del PE2[ID]
                else:
                    PE1[ID] = r1
            elif (flag & 0x80):  # If read 2 (last segment in template)
                if (flag & 0x10):  # if RevComp
                    r2 = ['R', line2[2], int(line2[3]) + cigar.getAlignmentLength() - 1, line2[3]]
                else:
                    r2 = ['F', line2[2], line2[3], int(line2[3]) + cigar.getAlignmentLength() - 1]
                if ID in PE1:
                    write_pair_tab(PE1[ID], r2)
                    del PE1[ID]
                else:
                    PE2[ID] = r2


print "Records processed: %s, discarded secondary alignments: %s, single end reads: %s, paired end reads: %s, broken pairs: %s" % (i / 2, secondary_alignment_count, single_end_count, paired_end_count, paired_end_broken)

# Finally go through PE1 and PE2, write out any SE reads that might exist:
# print "Checking for unmatched, paired reads and writing as SE"
# print "PE1 reads: ", len(PE1)
# print "PE2 reads: ", len(PE2)

# for k in PE1.keys():
#     outSE.write("@" + k + '#_1\n')
#     outSE.write(PE1[k][0] + '\n')
#     outSE.write('+\n' + PE1[k][1] + '\n')
# for k in PE2.keys():
#     outSE.write("@" + k + '#_2\n')
#     outSE.write(PE2[k][0] + '\n')
#     outSE.write('+\n' + PE2[k][1] + '\n')

insam.close()
output.close()
