#!/usr/bin/env python

import os, sys

ROOTimported = False

try:
    import ROOT
    ROOTimported = True

except ImportError:
    if os.environ.has_key('ROOTSYS'):
        sys.path.append(os.path.join(os.environ['ROOTSYS'],'lib'))

        try:
            import ROOT
            ROOTimported = True
        except ImportError:
            pass
        

if not ROOTimported:
    print >> sys.stderr,"could not import ROOT (pyROOT) module"
    print >> sys.stderr,"you could e.g. try to initialize a CMSSW"
    print >> sys.stderr,"runtime environment (using the command 'cmsenv')"
    print >> sys.stderr,"before starting this script."
    sys.exit(1)

#----------------------------------------------------------------------
def findConsecutiveRanges(values):

    retval = []

    if not values:
        return retval

    # at least one element
    values = sorted(values)

    lastValue = None

    # start a new range
    currentRange = [ ]

    for value in values:

        if not currentRange:
            currentRange.append(value)
            lastValue = value
            continue

        if value == lastValue + 1:
            lastValue = value
            continue


        # non-consecutive, close previous range
        currentRange.append(lastValue)

        retval.append(currentRange)
        currentRange = [ value ]

        lastValue = value

        continue

    assert(len(currentRange) != 2)

    # check whether we have started a range
    if len(currentRange) == 1:
        currentRange.append(lastValue)
        retval.append(currentRange)
        currentRange = []


    return retval



#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------    

ARGV = sys.argv[1:]

# first index is run number, second index is the lumi section,
# value is 1 (we assume that we're still running with python2.4
# which does not have sets...)
lumiSectionsFound = {}

for index, fname in enumerate(ARGV):
    print >> sys.stderr,"opening file",(index+1),"of",len(ARGV)
    fin = ROOT.TFile.Open(fname)

    event = fin.Get("event")

    event.SetEstimate(event.GetEntries())

    event.Draw("run:lumis","","goff")

    num = event.GetSelectedRows()
    runs = event.GetV1()
    lumis = event.GetV2()

    for i in range(num):
        run = int(runs[i] + 0.5)
        lumi = int(lumis[i] + 0.5)

        # print "run=",run,"lumi=",lumi
        lumiSectionsFound.setdefault(run,{})[lumi]=1

    fin.Close()


# from pprint import pprint
# pprint(lumiSectionsFound)

#----------------------------------------
# print what we've found
#----------------------------------------

# print findConsecutiveRanges([1,2,5,6,7,8,11,12])
# print findConsecutiveRanges([12])
# sys.exit(1)

runLines = []
for run in sorted(lumiSectionsFound.keys()):
    line = '"%d": ' % run

    line += str(findConsecutiveRanges(lumiSectionsFound[run].keys()))

    runLines.append(line)
print "{" + ", ".join(runLines) + "}"
