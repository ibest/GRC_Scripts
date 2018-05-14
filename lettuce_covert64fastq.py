#!/usr/bin/env python

"""
Very simply script which compares two fastq files

# Copyright 2015, Matt Settles
# Modified Aug 15, 2015
"""


import gzip
from Bio import SeqIO
import sys


fq1 = '101005_I125_FC808DGABXX_L1_LETabuDADDBBAPEI-12_1.fq.gz'
fq2 = '101005_I125_FC808DGABXX_L1_LETabuDADDBBAPEI-12_2.fq.gz'

fq1_out = 'Lsat-Davis_R1_001.fastq.gz'
fq2_out = 'Lsat-Davis_R2_001.fastq.gz'


if fq1.split('.')[-1] == 'gz':
    h1_in = SeqIO.parse(gzip.open(fq1, 'rb'), 'fastq-illumina')
else:
    h1_in = SeqIO.parse(open(fq1, 'r'), 'fastq-illumina')

if fq2.split('.')[-1] == 'gz':
    h2_in = SeqIO.parse(gzip.open(fq2, 'rb'), 'fastq-illumina')
else:
    h2_in = SeqIO.parse(open(fq2, 'r'), 'fastq-illumina')



if fq1_out.split('.')[-1] == 'gz':
    h1_out = gzip.open(fq1_out, 'ab')
else:
    h1_out = open(fq1_out, 'a')

if fq2_out.split('.')[-1] == 'gz':
    h2_out = gzip.open(fq2_out, 'ab')
else:
    h2_out = open(fq2_out, 'a')

#SEQ @B808DGABXX:1:1:1324:2041#NNTGTAAT/1
#MOD @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG

i = 0
same = 0
different = 0

fixed_r1 = []
fixed_r2 = []

try:
    while True:
        i += 1
        if i % 100000 == 0 and i > 0:
            print "records: %s \n\tidentical records: %s  \n\tdifferent records: %s" % (i, same, different)

        r1 = h1_in.next()
        r2 = h2_in.next()
        r1_id = r1.id.split("#")
        r2_id = r2.id.split("#")
        if str(r1_id[0]) != str(r2_id[0]):
            print "Mismatch record %s:" % i
            print '\t' + str(r1_id[0])
            print '\t' + str(r2_id[0])
            different += 1
        else:
            same += 1
        new_id = "UNK:1:" + r1_id[0][1:]
        r1.id = new_id
        r2.id = new_id
        r1.description = '1:N:0:'+ r1_id[1].split('/')[0]
        r2.description = '2:N:0:'+ r2_id[1].split('/')[0]
        
        fixed_r1.append(r1)
        fixed_r2.append(r2)
        
except StopIteration:
    pass
finally:
    SeqIO.write(fixed_r1,fq1_out,"fastq")
    SeqIO.write(fixed_r2,fq2_out,"fastq")
    #Print final output
    print "Finished processing one pair of files."
    print "Total records: %s \n\tidentical records: %s  \n\tdifferent records: %s" % (i, same, different)
