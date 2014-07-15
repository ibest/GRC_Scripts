#!/usr/bin/env python

"""
A script that checks whether a run is finished.
Inputs:
    RunFolder
outputs:
    if run is finished, create empty file "runfinished" in RunFolder

Overview:
    1. Read through RunInfo.xml, get the number of total cycles.
    2. Check that RunFolder/Data/Intensities/L001/CNNN.1

Alternative:
    1) To check whether your run has completed successfully or failed, you can check the end of the AnalysisLog.txt file in the RunInfo folder.
    The output would look similar to this:

    6/15/2012,10:18:28.384,Saving Completed Job Information
    6/15/2012,10:18:28.439,Copying Remaining Files To Network
    6/15/2012,10:21:04.397,Ending Execution for Analysis D:\Illumina\MiSeqAnalysis\110617_M3_0307_AFCA012C_1_Amplicon
    6/15/2012,10:21:04.413,Worker Thread Exit for Analysis D:\Illumina\MiSeqAnalysis\110617_M3_0307_AFCA012C_1_Amplicon
"""
import os
import sys

####### Functions #############
def pxml(l):
    kvp = {}
    l2 = l.strip().strip('<').strip('>')
    l3 = l2.split()
    for kv in l3:
        kv = kv.split('=')
        if len(kv) == 2:
            kvp[kv[0]] = kv[1].strip('"')
    return(kvp)
#######################################################


#First, get command line arguments:
if len(sys.argv) < 2:
    print "Error, could not find input parameter"
    print "Usage:"
    print "\tcheck_status.py RunFolder"
    sys.exit()

runfolder = os.path.realpath(sys.argv[1])

#Check that this is a real path:
if not os.path.exists(runfolder):
    print "Error, could not find RunFolder: %s\nExiting...." % runfolder
    sys.exit()

### process RunInfo.xml for number of cycles:  ###
infn = os.path.join(runfolder, "RunInfo.xml")
if os.path.isfile(infn):
    inf1 = open(infn, 'r')
else:
    print("Error cannot find %s" % inf1)
    sys.exit()

NumCycles = 0
TileCount = ''

for l in inf1:
    l = l.strip()
    if 'FlowcellLayout' in l:
        m = pxml(l)
        TileCount = m['TileCount']
    if 'NumCycles' in l:
        m = pxml(l)
        NumCycles += int(m['NumCycles'])
inf1.close()

print "Run has %s cycles and %s tiles." % (NumCycles, TileCount)

#Check whether the final file exists for NumCycles:
finalf = os.path.join(runfolder, "Data/Intensities/BaseCalls/L001/", ("C" + str(NumCycles) + ".1"), ("s_1_11" + TileCount.zfill(2) + ".bcl"))
finalf = os.path.realpath(finalf)
if os.path.isfile(finalf):
    print "Final file %s exists." % finalf
    out = open(os.path.join(runfolder, "runfinished"), 'w')
    out.write(" ")
    out.close()
else:
    print "Could not find file: %s." % finalf
