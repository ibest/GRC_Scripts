#!/usr/bin/env python

"""
Demultiplex reads when they come with a "barcode" read and haven't been split previously.
"""

import sys
from Bio import SeqIO
import os
import time


#------------------- classes -------------------------------

class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle

    def __iter__(self):
        return self

    def next(self):
        with self.inf as inf:
            lines = {'id': inf.readline().strip(),
                     'seq': inf.readline().strip(),
                     '+': inf.readline().strip(),
                     'qual': inf.readline().strip()}
            assert(len(lines['seq']) == len(lines['qual']))
        if lines['id'] == '' or lines['seq'] == '' or lines['+'] == '' or lines['qual'] == '':
            return False
        else:
            return lines

    @staticmethod
    def parse(handle):
        return fastqIter(handle)

    def close(self):
        self.inf.close()

"""
def writePE(outf1, outf2, r1, r2):
    #read1
    outf1.write(r1['id'])
    outf1.write(r1['seq'])
    outf1.write(r1['+'])
    outf1.write(r1['qual'])

    #read2
    outf2.write(r2['id'])
    outf2.write(r2['seq'])
    outf2.write(r2['+'])
    outf2.write(r2['qual'])

    outPE2.write("@" + ID + "#0/2" '\n')
    outPE2.write(r2[0] + '\n')
    outPE2.write('+\n' + r2[1] + '\n')
"""


#------------------- functions ------------------------------
def strdist(s1, s2):
    if len(s1) == len(s2) and len(s1) > 0:
        return sum(map(lambda x: x[0] != x[1], zip(s1, s2)))  # python is bad-ass
    else:
        print "ERROR lengths of barcodes and index read do not match!"
        print "Target", s1
        print "Index read:", s2
        sys.exit()

#-------------------- Testing --------------------------------

## Testing stuff:
# sys.argv.append("targets.tsv")
# sys.argv.append("DoubleBarcodeRun_140407_NoIndex_L001_R1_001.fastq.gz")
# sys.argv.append("DoubleBarcodeRun_140407_NoIndex_L001_R4_001.fastq.gz")
# sys.argv.append("DoubleBarcodeRun_140407_NoIndex_L001_R2_001.fastq.gz")
# sys.argv.append("DoubleBarcodeRun_140407_NoIndex_L001_R3_001.fastq.gz")


#----------- introductory stuff to check input-----------------
if len(sys.argv) not in [4, 5, 6]:
    print "Usage:"
    print "\tFor Dual index PE: Illumina_demultiplexer.py targets.tsv read1.fastq<.gz> read2.fastq<.gz> barcode1.fastq<.gz> barcode2.fastq<.gz>"
    print "\tFor PE: Illumina_demultiplexer.py targets.tsv read1.fastq<.gz> read2.fastq<.gz> barcode.fastq<.gz>"
    print "\tFor SE: Illumina_demultiplexer.py targets.tsv read1.fastq<.gz> barcode.fastq<.gz>"
    print "\t Where: "
    print "\t\tfastq files can be gzipped or not but must end in .gz if they are gzipped."
    print "\t\ttargets.tsv must be TAB seperated, and have column names sample_id<tab>barcode"
    print "\t\tor for dual index: sample_id<tab>barcode1<tab>barcode2"
    sys.exit()

#--------Get file names--------------
dualbc = False
PE = False

targetsf = sys.argv[1]
PE1f = sys.argv[2]

if len(sys.argv) == 4:
    barcodef1 = sys.argv[3]
elif len(sys.argv) == 5:
    PE2f = sys.argv[3]
    barcodef1 = sys.argv[4]
    PE = True
elif len(sys.argv) == 6:
    PE2f = sys.argv[3]
    barcodef1 = sys.argv[4]
    barcodef2 = sys.argv[5]
    dualbc = True
    PE = True

#--------- Check that all inputs exist, report error and exit if not -----------
if not os.path.isfile(targetsf):
    print targetsf, " file cannot be found, please check input parameters and try again."
    sys.exit()


if not os.path.isfile(PE1f):
    print PE1f, " file cannot be found, please check input parameters and try again."
    sys.exit()

if not os.path.isfile(barcodef1):
    print barcodef1, " file cannot be found, please check input parameters and try again."
    sys.exit()

if dualbc and not os.path.isfile(barcodef2):
    print barcodef2, " file cannot be found, please check input parameters and try again."
    sys.exit()

if PE and not os.path.isfile(PE2f):
    print PE2f, " file cannot be found, please check input parameters and try again."
    sys.exit()

#-----------Create file handles for Fastq/gz files -------------------
if PE1f.split(".")[-1] == "gz":
    print "Input appears to be gzipped."
    import gzip
    iter1 = SeqIO.parse(gzip.open(PE1f, 'rb'), "fastq")
    if PE:
        iter2 = SeqIO.parse(gzip.open(PE2f, 'rb'), "fastq")
    iteridx1 = SeqIO.parse(gzip.open(barcodef1, 'rb'), "fastq")
    if dualbc:
        iteridx2 = SeqIO.parse(gzip.open(barcodef2, 'rb'), "fastq")
else:
    print "Files appear to be uncompressed fastq."
    iter1 = SeqIO.parse(open(PE1f, 'r'), "fastq")
    if PE:
        iter2 = SeqIO.parse(open(PE2f, 'r'), "fastq")
    iteridx = SeqIO.parse(open(barcodef1, 'r'), "fastq")
    if dualbc:
        iteridx2 = SeqIO.parse(barcodef2, "fastaq")

#---------Open and parse the targets file -------------
targets = {}  # barcode is key, sample_id is value

with open(targetsf, 'r') as inf:
    l = inf.next()
    l = l.strip().split('\t')

    if dualbc:
        if l[0] != 'sample_id' or l[1] != 'barcode1' or l[2] != 'barcode2':
            print "Targets file is not properly formatted, first line must be: sample_id<tab>barcode1<tab>barcode2<additional field optional>"
            sys.exit()
    else:
        if l[0] != 'sample_id' or l[1] != 'barcode':
            print "Targets file is not properly formatted, first line must be: sample_id<tab>barcode<additional field optional>"
            sys.exit()

    for l in inf:
        if l == '\n':
            break
        l = l.strip().split("\t")
        if dualbc:
            k = l[1] + l[2]
        else:
            k = l[1]
        targets[k] = l[0]

# ------- Open output files -----------
outfiles = {}

for k in targets:
    if PE:
        outfiles[k] = [gzip.open("%s_%s_R1.fastq.gz" % (targets[k], k), 'wb'), gzip.open("%s_%s_R2.fastq.gz" % (targets[k], k), 'wb')]
    else:
        outfiles[k] = [gzip.open("%s_%s_R2.fastq.gz" % (targets[k], k), 'wb')]


# ------- make some counters for reporting ------------
counters = {}
for k in targets:
    counters[k] = [0, 0]  # first value is perfect matches, second is 1bp mismatch

counters['undetermined'] = 0

reads = 0
# Now iterate through both input files, split and write out
try:
    startt = time.time()
    while 1:
        read1 = iter1.next()
        if PE:
            read2 = iter2.next()
        idx1 = iteridx1.next()
        if dualbc:
            idx2 = iteridx2.next()

        reads += 1
        key = None
        k = idx1.seq.tostring()
        if dualbc:
            k += idx2.seq.tostring()

        if k in outfiles:
            key = k
            counters[key][0] += 1
        else:  # handle 1 base mismatches gracefully
            for bc in outfiles:
                if strdist(bc, k) <= 1:
                    key = bc
                    counters[key][1] += 1
        if key is not None:
            read1.description += key
            SeqIO.write(read1, outfiles[key][0], "fastq")
            if PE:
                read2.description += key
                SeqIO.write(read2, outfiles[key][1], "fastq")
        else:
            counters['undetermined'] += 1

        if reads % 10000 == 0:
            print"-------------------"
            print "%s reads processed in %s seconds, %s reads/sec" % (reads, (time.time() - startt), reads/(time.time() - startt))
            for k in targets:
                print "%s_%s: perfect match:%s, mismatch:%s" % (targets[k], k, counters[k][0], counters[k][1])
            print "undetermined: %s" % counters['undetermined']
            print"-------------------"

except StopIteration:
    pass
finally:
    print "Finished processing "
    print"-------------------"
    print "%s reads processed."
    for k in targets:
        print "%s_%s: perfect match:%s, mismatch:%s" % (targets[k], k, counters[k][0], counters[k][1])
    print "undetermind: %s" % counters['undetermined']
    print"-------------------"
