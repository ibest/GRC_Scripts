#!/usr/bin/env python
"""
# Copyright 2013, Sam Hunter, Brice Sarver, Matt Settles
# Modified Aug 27, 2014
"""


#from Bio import SeqIO
from optparse import OptionParser
import sys, os, os.path, time, gzip
from collections import Counter
from subprocess import Popen, PIPE, STDOUT
import string


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


class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle

    def __iter__(self):
        return self

    def next(self):
        lines = {'id': self.inf.readline().strip()[1:],
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
    handle.write('@' + fq['id'] + '\n')
    handle.write(fq['seq'] + '\n')
    handle.write(fq['+'] + '\n')
    handle.write(fq['qual'] + '\n')


## Parse options and setup ##
usage = "usage %prog -d [path to directory of raw reads] -o [output file prefix (path + name)]"
usage += "\n\te.g. %prog read1.fastq read2.fastq ./out/output_base"
parser = OptionParser(usage=usage)

parser.add_option('-d', '--directory', help="Directory containing read files to de-duplicate",
                  action="store", type="str", dest="sample_dir")

parser.add_option('-o', '--output', help="Directory + prefix to output de-duplicated reads",
                  action="store", type="str", dest="output_dir")

parser.add_option('-s', '--skip_dup', help="Skip de-dupping, merge files only and format for further processing in seqyclean",
                  action="store_true", dest="skip",default=False)


(options, args) = parser.parse_args()

sample_dir = options.sample_dir
output_dir = options.output_dir
skip = options.skip

if len(args) != 0 or sample_dir is None or output_dir is None:
    parser.print_help()
    sys.exit()


#kindly provided by http://stackoverflow.com/questions/7099290/how-to-ignore-hidden-files-using-os-listdir-python
#glob.glob will list hidden files
#this replaces that functionality when hidden files exist, like in my reads from Berkeley
def listdir_nohidden(path):
    for f in sorted(os.listdir(path), key=str.lower):
        if not f.startswith('.'):
            yield f


def main(infile1, infile2, outfile1, outfile2,skip):
#Set up the global variables
    global count
    global i
    global duplicates
    global rev
    global start
#Open inputs:
    if infile1.split(".")[-1] == "gz":
        #import gzip
        #iterator1 = SeqIO.parse(gzip.open(infile1, 'rb'), 'fastq')
        #iterator2 = SeqIO.parse(gzip.open(infile2, 'rb'), 'fastq')
        iterator1 = fastqIter(sp_gzip_read(infile1))
        iterator2 = fastqIter(sp_gzip_read(infile2))
    else:
        iterator1 = fastqIter(infile1)
        iterator2 = fastqIter(infile2)

    # elif infile1.split(".")[-1] == "fastq":
    #     #iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
    #     #iterator2 = SeqIO.parse(open(infile2, 'r'), 'fastq')
    #     iterator1 = fastqIter(infile1)
    #     iterator2 = fastqIter(infile2)

    # else:
    #     iterator1 = SeqIO.parse(open(infile1, 'r'), 'fastq')
    #     iterator2 = SeqIO.parse(open(infile2, 'r'), 'fastq')

    try:
        while 1:
            c = 0
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            if skip is True: ## skip dedup, just write out
                writeFastq(outfile1, seq1)
                writeFastq(outfile2, seq2)
            else:
                comb = seq1['seq'][10:35] + seq2['seq'][10:35]
                rcomb = revcomp(comb)
                if rcomb in count:
                    count[rcomb] += 1
                    c = count[rcomb]
                    rev += 1
                else:
                    count[comb] += 1
                    c = count[comb]
                if c == 1:
                    writeFastq(outfile1, seq1)
                    writeFastq(outfile2, seq2)
                else:
                    duplicates += 1
            i += 1
            if i % 100000 == 0:
                print "Pairs:", "| reads:", i, "| duplicates:", duplicates, "| fwd:", duplicates-rev, "| rev:", rev, "| percent:", round(100.0*duplicates/i, 2), "| reads/sec:", round(i/(time.time() - start), 0)

    except StopIteration:
        pass
    finally:
#Print final output
        print "Finished processing one pair of files."


#####################################
count = Counter()
i = 0
duplicates = 0
rev = 0
start = time.time()

outfile1 = sp_gzip_write(output_dir + "_nodup_PE1.fastq.gz")
outfile2 = sp_gzip_write(output_dir + "_nodup_PE2.fastq.gz")

#outfile1 = gzip.open(output_dir + "_nodup_PE1.fastq.gz", 'wb')
#outfile2 = gzip.open(output_dir + "_nodup_PE2.fastq.gz", 'wb')

files = listdir_nohidden('./' + sample_dir)

print "skip detection of duplicates: " + str(skip)
for f in files:
    if "_R1" in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "_R2".join(f.split("_R1"))))
        main(infile1, infile2, outfile1, outfile2, skip)
    elif "READ1" in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "READ2".join(f.split("READ1"))))
        main(infile1, infile2, outfile1, outfile2, skip)
    elif "_1.fastq" in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "_2.fastq".join(f.split("_1.fastq"))))
        main(infile1, infile2, outfile1, outfile2, skip)
    elif '_PE1.fastq' in f:
        print f
        infile1 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, f))
        infile2 = os.path.realpath(os.path.join(os.getcwd(), sample_dir, "_PE2.fastq".join(f.split("_PE1.fastq"))))
        main(infile1, infile2, outfile1, outfile2, skip)
    else:
        print "%s not recognized" % f

outfile1.close()
outfile2.close()

print "Final:", "| reads:", i, "| duplicates:", duplicates, "| fwd:", duplicates-rev, "| rev:", rev, "| percent:", round(100.0*duplicates/i, 2), "| reads/sec:", round(i/(time.time() - start), 0)
