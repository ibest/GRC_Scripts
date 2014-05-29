#!/usr/bin/env python

"""
Parse output from FLASH
"""

###
#[FLASH] Parameters:
#[FLASH]     Min overlap:          10
#[FLASH]     Max overlap:          65
#[FLASH]     Phred offset:         33
#[FLASH]     Combiner threads:     40
#[FLASH]     Max mismatch density: 0.250000
#[FLASH]     Output format:        gzip
#[FLASH]     Interleaved input:    false
#[FLASH]     Interleaved output:   false
#[FLASH]
#[FLASH] Starting FASTQ readers and writer threads
#                   ...
#[FLASH]
#[FLASH] Read combination statistics:
#[FLASH]     Total reads:      8004022
#[FLASH]     Combined reads:   4880119
#[FLASH]     Uncombined reads: 3123903
#[FLASH]     Percent combined: 60.97%
###

import sys
import os
import re


def parse_flash(fileinput_stream):
    skip = 4
    result = {}
    for i, line in enumerate(fileinput_stream):
        if skip == 4:
            ### parse version
            result['Flash_version'] = re.split(r' +',line.rstrip())[3]
            skip = 3
            continue
        if skip == 3 and not "Parameters" in line:
            continue
        elif skip == 3 and "Parameters" in line:
            skip = 2
            continue
        elif skip == 2 and not "Starting FASTQ readers and writer threads" in line:
            ### parse Parameters
            data = re.split(': +', re.sub(r'\[FLASH\] +','',line.rstrip()))
            if len(data) == 2:
                name = re.sub(r' ','_',data[0])
                result[name] = data[1]
            continue
        elif skip == 2 and "Starting FASTQ readers and writer threads" in line:
            skip = 1
            continue
        elif skip == 1 and not "Read combination statistics" in line:
            continue
        elif skip == 1 and  "Read combination statistics" in line:
            skip =  0
            continue
        elif skip == 0 and not "Writing histogram files" in line:
            ### parse read combination statistics
            data = re.split(': +', re.sub(r'\[FLASH\] +','',line.rstrip()))
            if len(data) == 2:
                name = re.sub(r' ','_',data[0])
                result[name] = data[1]
            continue
        elif skip == 0 and "Writing histogram files" in line:
            return(result)

if __name__ == '__main__':
    # grab command line args if any
    if (len(sys.argv) > 2):
        sys.exit('Usage: %s [flash output file]' % sys.argv[0])
    elif (len(sys.argv) == 2):
        # file provided on command line
        if not os.path.exists(sys.argv[1]):
            sys.exit('ERROR: file %s was not found!' % sys.argv[1])
        with open(sys.argv[1], 'rb') as infile:
            output = parse_flash(infile)
            print output
    else:
        # streaming via standard in
        with sys.stdin as infile:
            output = parse_flash(infile)
            print output
