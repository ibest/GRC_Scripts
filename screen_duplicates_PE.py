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


## Parse options and setup ##
usage = "usage %prog -d [path to directory of raw reads] -o [output file prefix (path + name)] -(blsa) --quite"
usage += "screen_duplicates_PE.py will process all pairs in the provided directory"
usage += "\n\tif no directory is provided read1 and read2 must be supplied on the command line"
usage += "\n\te.g. %prog read1.fastq read2.fastq"
parser = OptionParser(usage=usage,version="%prog 2.0.0")

parser.add_option('-d', '--directory', help="Directory containing read files to de-duplicate",
                  action="store", type="str", dest="sample_dir", default=None)

parser.add_option('-o', '--output', help="Directory + prefix to output de-duplicated reads",
                  action="store", type="str", dest="output_dir", default="reads")

parser.add_option('-b', '--start', help="position to start duplication check",
                  action="store", type="int", dest="start", default=10)

parser.add_option('-l', '--length', help="length of duplication check",
                  action="store", type="int", dest="length", default=25)

parser.add_option('-s', '--skip_dup', help="Skip de-dupping, merge files only and format for further processing in seqyclean",
                  action="store_true", dest="skip",default=False)

parser.add_option('-a', '--sra', help="Data was downloaded from the SRA, requires ID's to be rewritten",
                  action="store_true", dest="sra",default=False)

parser.add_option('--quite', help="turn off verbose output",
                  action="store_false", dest="verbose",default=True)


(options, args) = parser.parse_args()

sample_dir = options.sample_dir
output_dir = options.output_dir
skip = options.skip

start = options.start - 1
end = start + options.length

#if len(args) != 0 or sample_dir is None or output_dir is None:
#    parser.print_help()
#    sys.exit()

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
    global stime
#Open inputs:
    if infile1.split(".")[-1] == "gz":
        iterator1 = fastqIter(sp_gzip_read(infile1))
        iterator2 = fastqIter(sp_gzip_read(infile2))
    else:
        iterator1 = fastqIter(infile1)
        iterator2 = fastqIter(infile2)
    try:
        while 1:
            c = 0
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            if skip is True: ## skip dedup, just write out
                writeFastq(outfile1, seq1)
                writeFastq(outfile2, seq2)
            else:
                comb = seq1['seq'][start:end] + seq2['seq'][start:end]
                count[comb] += 1
                c = count[comb]
                if c == 1:
                    if options.sra: ## modify read id adding in index of read and pair information needed
                        seq1['id'] = "@HWI-"+i+":0:0:0:0:0:0 1:Y:0:"
                        seq2['id'] = "@HWI-"+i+":0:0:0:0:0:0 2:Y:0:"
                    writeFastq(outfile1, seq1)
                    writeFastq(outfile2, seq2)
                else:
                    duplicates += 1
            i += 1
            if i % 100000 == 0 and options.verbose:
                print "Pairs:", "| reads:", i, "| duplicates:", duplicates, "| percent:", round(100.0*duplicates/i, 2), "| reads/sec:", round(i/(time.time() - stime), 0)

    except StopIteration:
        pass
    finally:
#Print final output
        print "Finished processing one pair of files."


#####################################
count = Counter()
i = 0
duplicates = 0
stime = time.time()

outfile1 = sp_gzip_write(output_dir + "_nodup_PE1.fastq.gz")
outfile2 = sp_gzip_write(output_dir + "_nodup_PE2.fastq.gz")

if sample_dir is not None:
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
else: ## files on the command line
    if len(args) != 2:
        parser.error("incorrect number of arguments, expecting 2 reads files")
    print "skip detection of duplicates: " + str(skip)
    infile1 = args[0]
    infile2 = args[1]
    main(infile1, infile2, outfile1, outfile2, skip)


outfile1.close()
outfile2.close()

print "Final:", "| reads:", i, "| duplicates:", duplicates, "| percent:", round(100.0*duplicates/i, 2), "| reads/sec:", round(i/(time.time() - stime), 0)
