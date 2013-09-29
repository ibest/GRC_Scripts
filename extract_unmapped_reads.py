#!/usr/bin/env python

'''
Iterate over a SAM or SAM.gz file, take everything where the 3rd and
4th flag bit are set to 1 and write reads out to files.
The SAM definition requires that pairs appear consequitively for PE reads,
so the trick is to read in two lines, compare them to make sure
they are paired, and then output to .fastq.

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
'''
import sys
import os
#from collections import OrderedDict
from optparse import OptionParser  # http://docs.python.org/library/optparse.html
import gzip


usage = "usage: %prog [options] inputfile.SAM output_base"
parser = OptionParser(usage=usage)
parser.add_option('-u', '--uncompressed', help="leave output files uncompressed",
                  action="store_true", dest="uncompressed")

(options,  args) = parser.parse_args()  # uncomment this line for command line support

if len(args) < 1:
    parser.print_help()
    sys.exit()

infile = args[0]
base = args[1]

PE1 = {}
PE2 = {}

#Start opening input/output files:
if not os.path.exists(infile):
    print "Error, can't find input file %s" % infile
    sys.exit()

if infile.split(".")[-1] == "gz":
    insam = gzip.open(infile, 'rb')
else:
    insam = open(infile, 'r')

if options.uncompressed:
    outPE1 = open(base + "_PE1.fastq", 'w')
    outPE2 = open(base + "_PE2.fastq", 'w')
    outSE = open(base + "_SE.fastq", 'w')
else:
    outPE1 = gzip.open(base + "_PE1.fastq.gz", 'wb')
    outPE2 = gzip.open(base + "_PE2.fastq.gz", 'wb')
    outSE = gzip.open(base + "_SE.fastq.gz", 'wb')


def writeread(ID, r1, r2):
    #read1
    outPE1.write("@" + ID + "#0/1" '\n')
    outPE1.write(r1[0] + '\n')
    outPE1.write('+\n' + r1[1] + '\n')
    #read2
    outPE2.write("@" + ID + "#0/2" '\n')
    outPE2.write(r2[0] + '\n')
    outPE2.write('+\n' + r2[1] + '\n')

i = 0
PE_written = 0
SE_written = 0
unwritten = 0
for line in insam:
    ## Check to see if it is a header line, if so skip and move on
    if line[0] == "@":
        continue
    ## not a header line, so increment reads seen
    line2 = line.strip().split()
    if len(line2) > 2:
        flag = int(line2[1])
    else:
        sys.exit("Something wrong with line %s, does not contain at least 2 columns" % i+1)
    ## check for single end unmapped read
    #logic:  !(0x1) = multiple segments in sequencing,   0x4 = segment unmapped,
    # which means (is single, and this segment is unmapped)
    if flag & 0x1 == 0 and flag & 0x4 != 0:
        ID = line2[0].split("#")[0]
        outSE.write("@" + ID + '\n')
        outSE.write(line2[9] + '\n')
        outSE.write('+\n' + line2[10] + '\n')
        SE_written += 1
    ## check for both paired end unmapped read
    #logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
    # which means (is paired, and this segment is unmapped and the next segment is unmapped)
    elif flag & 0x1 != 0 and flag & 0x4 != 0 and flag & 0x8 != 0:
        if flag & 0x40 != 0:  # is this PE1 (first segment in template)
            #PE1 read, check that PE2 is in dict and write out
            ID = line2[0].split("#")[0]
            r1 = [line2[9], line2[10]]
            if ID in PE2:
                writeread(ID, r1, PE2[ID])
                del PE2[ID]
                PE_written += 1
            else:
                PE1[ID] = r1
        elif flag & 0x80 != 0:  # is this PE2 (last segment in template)
            #PE2 read, check that PE1 is in dict and write out
            ID = line2[0].split("#")[0]
            r2 = [line2[9], line2[10]]
            if ID in PE1:
                writeread(ID, PE1[ID], r2)
                del PE1[ID]
                PE_written += 1
            else:
                PE2[ID] = r2
    else:
        unwritten += 1
    i += 1
    if i % 100000 == 0:
        print "Record %s" % i
        print "\t Unwritten: %s" % unwritten
        print "\t Written: PE:%s SE:%s" % (PE_written,SE_written)

## Finally go through PE1 and PE2, write out any SE reads that might exist:
print "Checking for unmatched, paired reads"
print "PE1 reads: ", len(PE1)
print "PE2 reads: ", len(PE2)

outPE1.close()
outPE2.close()
outSE.close()