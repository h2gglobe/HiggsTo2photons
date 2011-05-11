from optparse import OptionParser
import os
from os import popen, path, environ, getpid
import re
import data_replica
from ROOT import *

parser = OptionParser()
parser.add_option("-c", "--clean", action="store_true", dest="clean", default=False, help="Remove Corrupted File")
parser.add_option("-d", "--duplicate", action="store_true", dest="duplicate", default=False, help="Remove Duplicate Job Submissions")
parser.add_option("-p", "--path", dest="directory", help="The data directory on CASTOR", metavar="DIR")
parser.add_option("-r", "--recursive", action="store_true", dest="recursive", default=False, help="Transfers files in sub directories")
parser.add_option("--single-thread", action="store_true", dest="singlethread", default=False, help="Disables parallel transfers for lcg-cp")
parser.add_option("--to-site", action="store", dest="TO_SITE", default="", help="Destination file, eg: T2_CH_CSCS")
parser.add_option("--remote-directory", action="store", dest="USERDIR", default="", help="Destination directory, eg: /store/user/drberry")
parser.add_option("--dryrun", action="store_true", dest="DRYRUN", default=False, help="Do not actually copy anything")
parser.add_option("--castor-stage", action="store_true", dest="CASTORSTAGE", default=False, help="Enables staging of Castor files in a local tmp dir. Works only on lxplus, and uses $TMPDIR as tmp dir.")
(options, args) = parser.parse_args()

### passing options to data_replica
class drOptions:
	usePDS = False
	Replicate = False
	RECREATE_SUBDIRS = False
	CASTORSTAGE = options.CASTORSTAGE
	DEBUG = False
	TOOL='lcg-cp'
	DRYRUN = options.DRYRUN
	pass

def cleanfiles(dir):
		filename = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()
		for i in range(len(filename)):
			filename[i] = filename[i].strip("\n")
			testfile = TFile.Open("rfio:"+dir+"/"+filename[i])
			if (testfile.IsZombie()):
				newfilename = filename[i].replace(".root",".resubmit")
				print "Moving corrupted file %s to %s" %(filename[i], jobnum[i], subnum[i])
				popen("rfrename "+dir+filename[j]+" "+dir+newfilename)
			else:
				TestTree = testfile.Get("event")
				if TestTree.GetEntries()==0:
					newfilename = filename[i].replace(".root",".empty")
					print "Moving corrupted file %s to %s" %(filename[i], jobnum[i], subnum[i])
					popen("rfrename "+dir+filename[j]+" "+dir+newfilename)

def removeduplicated(dir):
	jobnum=[]
	subnum=[]
	if dir[len(dir)-1:len(dir)]!="/": dir+="/" 
	filename = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()
	for i in range(len(filename)):
		filename[i] = filename[i].strip("\n")
		regsearch = re.search('[0-9]*_[0-9]*_[A-Za-z0-9]*.root',filename[i])

		if (regsearch):
			numarg = regsearch.group(0)
			numarg = numarg[0:numarg.rfind("_")]
			jobnum.append(int(numarg[0:numarg.find("_")]))
			subnum.append(int(numarg[numarg.find("_")+1:len(numarg)]))

	for i in range(max(jobnum)):
		if (jobnum.count(i)>1):
			lastsubmission=max(subnum[jobnum.index(i):jobnum.index(i)+jobnum.count(i)-1])
			firstsubmission=min(subnum[jobnum.index(i):jobnum.index(i)+jobnum.count(i)-1])
			if lastsubmission!=firstsubmission:
				for j in range(jobnum.index(i),jobnum.index(i)+jobnum.count(i)-1):
					newfilename = filename[j].replace(".root",".duplicate")
					print "Moving duplicate file %s to %s" %(filename[j],newfilename)
					popen("rfrename "+dir+filename[j]+" "+dir+newfilename)
			elif lastsubmission==firstsubmission:
				timesort = popen("rfdir "+dir+" | grep "+filename[jobnum.index(i)][:filename[jobnum.index(i)].rfind("_")]+" | sort -k 7,8 | awk '{print $9}' ").readlines()
				for j in range(len(timesort)-1):
					timesort[j] = timesort[j].strip("\n")
					newfilename = timesort[j].replace(".root",".duplicate")
					print "Moving duplicate file %s to %s" %(timesort[j],newfilename)
					popen("rfrename "+dir+timesort[j]+" "+dir+newfilename)

def makefilelists(dir):
	dirlist = dir.split("/")
	outputfile = dirlist[len(dirlist)-2]+".filelist"
	print "Making %s for directory %s" %(outputfile,dir)
	popen("nsfind "+dir+" | grep '.root\|.json' >& "+outputfile)

dir = options.directory
if (dir[len(dir)-1:len(dir)]!="/"): dir+="/" 
rootfiles = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()

if len(rootfiles)!=0:
	if (options.clean): cleanfiles(dir)
	if (options.duplicate): removeduplicated(dir)
	makefilelists(dir)

if options.recursive:
	subdirs = popen("rfdir "+dir+" | awk '{print $9}'").readlines()
	for i in range(len(subdirs)):
		subdir = dir+subdirs[i].strip("\n") + "/"
		if (options.clean): cleanfiles(subdir)
		if (options.duplicate): removeduplicated(subdir)
		makefilelists(subdir)

# Transfer Files
myOptions = drOptions()
myOptions.TO_SITE = options.TO_SITE
myOptions.FROM_SITE = "CERN_CASTOR_USER"

filelists = popen("ls | grep '\.filelist'").readlines()
for i in range(len(filelists)):
	filelist = filelists[i].strip("\n")
	myOptions.logfile = filelist[0:filelist.find(".")]+".log"
	fileName = [filelist,options.USERDIR+"/"+filelist[0:filelist.find(".")]]
	#print "fileName is: %s" %fileName
	NumFiles = popen("cat "+filelist+" | wc | awk '{print $1}'").readlines()
	NumFiles[0] = NumFiles[0].strip("\n")
	myOptions.THREADS = NumFiles[0]
	if int(NumFiles[0])>25: myOptions.THREADS="25"
	if options.singlethread: myOptions.THREADS="1"
	drExit = data_replica.data_replica(fileName, myOptions)
	if (drExit!=0 and not(options.DRYRUN)):
		print "Some errors in copying"
