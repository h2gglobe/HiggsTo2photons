from optparse import OptionParser
from os import popen
import re
import ROOT

parser = OptionParser()
parser.add_option("-c", "--clean", action="store_true", dest="clean", default=False, help="Remove Corrupted File")
parser.add_option("-k", "--check", action="store_true", dest="check", default=True, help="Check that all jobs are completed and that there are no duplicates.")
parser.add_option("-d", "--duplicate", action="store_true", dest="duplicate", default=False, help="Remove Duplicate Job Submissions")
parser.add_option("-p", "--path", dest="directory", help="The data directory on CASTOR", metavar="DIR")
parser.add_option("-r", "--recursive", action="store_true", dest="recursive", default=False, help="Clean files in sub directories")
parser.add_option("--eos", action="store_true", dest="eos", default=False, help="Enable EOS flag")
parser.add_option("--castor", action="store_true", dest="castor", default=False, help="Enable CASTOR flag")
parser.add_option("--debug", action="store_true", dest="debug", default=False, help="Enable Debug Mode")
parser.add_option("--dryrun", action="store_true", dest="dryrun", default=False, help="Dont Actually Move Files")
(options, args) = parser.parse_args()

if not options.eos and not options.castor:
    if options.directory.find("/castor/cern.ch/user/")!=-1: options.castor=True
    if options.directory.find("/eos/cms/store/")!=-1: options.eos=True

if (not options.eos and not options.castor) or (options.eos and options.castor):
    print "You must select either EOS xor CASTOR"
    exit(1)

eos="/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select"

def makejobstring(list):
    returnstring=""
    for i in list: returnstring+=str(i)+","
    return returnstring

def cleanfiles(dir):
    if dir[-1]!="/": dir+="/"
    print "Cleaning Directory:",dir
    if options.castor: popen("nsfind "+dir+" -type f | grep .root | xargs -i stager_qry -M {} | grep 'not on disk' | awk '{print $8}' | xargs -i stager_get -M {}")
    if options.castor: filename = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()
    if options.eos: filename = popen(eos+" ls "+dir+" | grep .root").readlines()
  
    if len(filename)==0: return
    for i in range(len(filename)):
        filename[i] = filename[i].strip("\n")
        if options.debug: print "Checking File: %s" %(dir+filename[i])
        if options.castor: testfile = ROOT.TFile.Open("rfio:"+dir+filename[i])
        if options.eos: testfile = ROOT.TFile.Open("root://eoscms/"+dir+filename[i])
        if (testfile==None or testfile.IsZombie()):
            if testfile!=None: testfile.Close()
            newfilename = filename[i].replace(".root",".resubmit")
            print "Moving corrupted file %s to %s" %(dir+filename[i], newfilename)
            if not options.dryrun:
                if options.castor: popen("rfrename "+dir+filename[i]+" "+dir+newfilename)
                if options.eos:
                    popen("xrdcp -s root://eoscms/"+dir+filename[i]+" root://eoscms/"+dir+newfilename)
                    popen(eos+" rm "+dir+filename[i])
        else:
            TestTree = testfile.Get("event")
            if TestTree.GetEntries()==0:
                testfile.Close()
                newfilename = filename[i].replace(".root",".empty")
                print "Moving empty file %s to %s" %(dir+filename[i], newfilename)
                if not options.dryrun:
                    if options.castor: popen("rfrename "+dir+filename[i]+" "+dir+newfilename)
                    if options.eos:
                        popen("xrdcp -s root://eoscms/"+dir+filename[i]+" root://eoscms/"+dir+newfilename)
                        popen(eos+" rm "+dir+filename[i])
            else: testfile.Close()

def checkfiles(dir):
    if dir[-1]!="/": dir+="/"
    print "Checking Directory:",dir
    if options.castor:
        print "Warning: Job checker not implemented for CASTOR!"
        return
    jobnumbers=[]
    missingjobs=[]
    duplicatejobs=[]
    filelist=popen(eos+" ls "+dir+" | grep .root").readlines()
    for file in filelist:
        jobnumber=int(file.strip("\n").split("_")[-3])
        jobnumbers.append(jobnumber)
    if len(jobnumbers)==0: return
    print max(jobnumbers),"total jobs found for",dir
    for jobnum in range(1,max(jobnumbers)+1):
        jobsearch=jobnumbers.count(jobnum)
        if jobsearch==0:
            print "Warning!!!! Job "+str(jobnum)+" not found!!!!"
            missingjobs.append(jobnum)
        if jobsearch>1:
            print "Warning!!!! Job "+str(jobnum)+" has "+str(jobsearch-1)+" duplicates!!!!"
            duplicatejobs.append(jobnum)
    if len(missingjobs)>0: print "Missing Jobs",makejobstring(missingjobs)
    if len(duplicatejobs)>0: print "Duplicate Jobs",makejobstring(duplicatejobs)

def removeduplicated(dir):
    if dir[-1]!="/": dir+="/"
    jobnum=[]
    subnum=[]
    print "Scanning for Duplicates in:",dir
    if options.castor: filename = popen("rfdir "+dir+" | grep .root | awk '{print $9}' ").readlines()
    if options.eos: filename = popen(eos+" ls "+dir+" | grep .root").readlines()
    if len(filename)==0: return
    for i in range(len(filename)):
        filename[i] = filename[i].strip("\n")
        regsearch = re.search('[0-9]*_[0-9]*_[A-Za-z0-9]*.root',filename[i])

        if regsearch:
            numarg = regsearch.group(0)
            numarg = numarg[0:numarg.rfind("_")]
            jobnum.append(int(numarg[0:numarg.find("_")]))
            subnum.append(int(numarg[numarg.find("_")+1:len(numarg)]))

    for i in range(max(jobnum)+1):
        if jobnum.count(i)>1:
            lastsubmission=max(subnum[jobnum.index(i):jobnum.index(i)+jobnum.count(i)])
            firstsubmission=min(subnum[jobnum.index(i):jobnum.index(i)+jobnum.count(i)])
            
            if lastsubmission!=firstsubmission:
                for j in range(jobnum.index(i),jobnum.index(i)+jobnum.count(i)-1):
                    newfilename = filename[j].replace(".root",".duplicate")
                    print "Moving duplicate file %s to %s" %(dir+filename[j],newfilename)
                    if not options.dryrun:
                        if options.castor: popen("rfrename "+dir+filename[j]+" "+dir+newfilename)
                        if options.eos:
                            popen("xrdcp -s root://eoscms/"+dir+filename[j]+" root://eoscms/"+dir+newfilename)
                            if options.debug: print eos+" rm "+dir+filename[j]
                            popen(eos+" rm "+dir+filename[j])
                            
            elif lastsubmission==firstsubmission:
                if options.castor: timesort = popen("rfdir "+dir+" | grep "+filename[jobnum.index(i)][:filename[jobnum.index(i)].rfind("_")]+" | sort -k 7,8 | awk '{print $9}' ").readlines()
                if options.eos: timesort = popen(eos+" ls -l "+dir+" | grep "+filename[jobnum.index(i)][:filename[jobnum.index(i)].rfind("_")]+" | sort -k 7,8 | awk '{print $9}' ").readlines()
                for j in range(len(timesort)-1):
                    timesort[j] = timesort[j].strip("\n")
                    newfilename = timesort[j].replace(".root",".duplicate")
                    print "Moving duplicate file %s to %s" %(dir+timesort[j],newfilename)
                    if not options.dryrun:
                        if options.castor: popen("rfrename "+dir+timesort[j]+" "+dir+newfilename)
                        if options.eos:
                            popen("xrdcp -s root://eoscms/"+dir+timesort[j]+" root://eoscms/"+dir+newfilename)
                            if options.debug: print eos+" rm "+dir+timesort[j]
                            popen(eos+" rm "+dir+timesort[j])

if options.debug: print "EOS Value:"+str(options.eos)
if options.debug: print "CASTOR Value:"+str(options.castor)

dir = options.directory
if dir[:1]!="/": dir+="/" 
if dir[0:1]!="/": dir="/"+dir 
if options.castor: rootfiles = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()
if options.eos: rootfiles = popen(eos+" ls "+dir+" | grep .root").readlines()

if len(rootfiles)!=0:
    if (options.duplicate): removeduplicated(dir)
    if (options.clean): cleanfiles(dir)
    if (options.check): checkfiles(dir)

if options.recursive:
    if options.castor: subdirs = popen("nsfind "+dir+" -type d ").readlines()
    if options.eos: subdirs = popen(eos+" find -d "+dir+" | grep /").readlines()
    for i in range(len(subdirs)):
        subdir = subdirs[i].strip("\n")
        if (options.duplicate): removeduplicated(subdir)
        if (options.clean): cleanfiles(subdir)
        if (options.check): checkfiles(subdir)
