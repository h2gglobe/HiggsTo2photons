from optparse import OptionParser
import os
from os import popen, path, environ, getpid
import re
from ROOT import *

parser = OptionParser()
parser.add_option("-c", "--clean", action="store_true", dest="clean", default=False, help="Remove Corrupted File")
parser.add_option("-d", "--duplicate", action="store_true", dest="duplicate", default=False, help="Remove Duplicate Job Submissions")
parser.add_option("-p", "--path", dest="directory", help="The data directory on CASTOR", metavar="DIR")
parser.add_option("-r", "--recursive", action="store_true", dest="recursive", default=False, help="Transfers files in sub directories")
(options, args) = parser.parse_args()

def cleanfiles(dir):
	if dir[len(dir)-1:len(dir)]!="/": dir+="/" 
	filename = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()
	for i in range(len(filename)):
		filename[i] = filename[i].strip("\n")
		testfile = TFile.Open("rfio:"+dir+filename[i])
		if testfile.IsZombie():
			newfilename = filename[i].replace(".root",".resubmit")
			print "Moving corrupted file %s to %s" %(filename[i], newfilename)
			popen("rfrename "+dir+filename[i]+" "+dir+newfilename)
		else:
			TestTree = testfile.Get("event")
			if TestTree.GetEntries()==0:
				newfilename = filename[i].replace(".root",".empty")
				print "Moving corrupted file %s to %s" %(filename[i], newfilename)
				popen("rfrename "+dir+filename[i]+" "+dir+newfilename)

def removeduplicated(dir):
	jobnum=[]
	subnum=[]
	if dir[len(dir)-1:len(dir)]!="/": dir+="/" 
	filename = popen("rfdir "+dir+" | grep .root | awk '{print $9}' ").readlines()
	for i in range(len(filename)):
		filename[i] = filename[i].strip("\n")
		regsearch = re.search('[0-9]*_[0-9]*_[A-Za-z0-9]*.root',filename[i])

		if regsearch:
			numarg = regsearch.group(0)
			numarg = numarg[0:numarg.rfind("_")]
			jobnum.append(int(numarg[0:numarg.find("_")]))
			subnum.append(int(numarg[numarg.find("_")+1:len(numarg)]))

	for i in range(max(jobnum)):
		if jobnum.count(i)>1:
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

dir = options.directory
if dir[len(dir)-1:len(dir)]!="/": dir+="/" 
rootfiles = popen("rfdir "+dir+" | awk '{print $9}' | grep .root").readlines()

if len(rootfiles)!=0:
	if (options.clean): cleanfiles(dir)
	if (options.duplicate): removeduplicated(dir)

if options.recursive:
	subdirs = popen("rfdir "+dir+" | awk '{print $9}'").readlines()
	for i in range(len(subdirs)):
		subdir = dir+subdirs[i].strip("\n") + "/"
		if (options.clean): cleanfiles(subdir)
		if (options.duplicate): removeduplicated(subdir)
