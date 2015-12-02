#!/usr/bin/env python

import sys
import os
import gzip


if len(sys.argv) == 2:
    infile = sys.argv[1]
    if not os.path.exists(infile):
        print "Error, can't find input file %s" % infile
        sys.exit()

    if infile.split(".")[-1] == "gz":
        insam = gzip.open(infile, 'rb')
    else:
        insam = open(infile, 'r')
else:
    ## reading from stdin
    insam = sys.stdin

output = os.path.join(os.path.dirname(os.path.abspath(fq1)),"softclip_inserts.log")

i = 0
PE_seen = 0
SE_seen = 0

for line in insam:
    if i % 100000 == 0 and i > 0 and options.verbose:
        print "Records processed: %s, PE_written: %s, SE_written: %s" % (i, PE_written, SE_written)
    #Comment/header lines start with @
    if line[0] != "@" and len(line.strip().split()) > 2:
        i += 1
        line2 = line.strip().split()
        flag = int(line2[1])
        #Handle SE:
        if (flag & 0x100): # secondary alignment
            continue
        # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
        if not (flag & 0x1) and not (flag & 0x4):
            ID = line2[0].split("#")[0]
            if (flag & 0x10):
                line2[9] = reverseComplement(line2[9])
                line2[10] = reverse(line2[10])
            outSE.write("@" + ID + '\n')
            outSE.write(line2[9] + '\n')
            outSE.write('+\n' + line2[10] + '\n')
            SE_written += 1
            continue
        #Handle PE:
        #logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
        if ((strict and (flag & 0x1) and not (flag & 0x4) and not (flag & 0x8))
                or (not strict and (flag & 0x1) and (not (flag & 0x4) or not (flag & 0x8)))):
            if (flag & 0x40):  # is this PE1 (first segment in template)
                #PE1 read, check that PE2 is in dict and write out
                ID = line2[0].split("#")[0]
                if (flag & 0x10):
                    line2[9] = reverseComplement(line2[9])
                    line2[10] = reverse(line2[10])
                r1 = [line2[9], line2[10]]  # sequence + qual
                if ID in PE2:
                    writeread(ID, r1, PE2[ID])
                    del PE2[ID]
                    PE_written += 1
                else:
                    PE1[ID] = r1
            elif (flag & 0x80):  # is this PE2 (last segment in template)
                #PE2 read, check that PE1 is in dict and write out
                ID = line2[0].split("#")[0]
                if (flag & 0x10):
                    line2[9] = reverseComplement(line2[9])
                    line2[10] = reverse(line2[10])
                r2 = [line2[9], line2[10]]
                if ID in PE1:
                    writeread(ID, PE1[ID], r2)
                    del PE1[ID]
                    PE_written += 1
                else:
                    PE2[ID] = r2

print "Records processed: %s, PE_written: %s, SE_written: %s" % (i, PE_written, SE_written)
## Finally go through PE1 and PE2, write out any SE reads that might exist:
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

outPE1.close()
outPE2.close()
outSE.close()
