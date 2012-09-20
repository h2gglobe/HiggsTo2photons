#!/usr/bin/python

import os
import sys
import pprint
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i","--input",dest="input_file",help="input txt file",default='temp_progress.txt')
parser.add_option("-b","--bash",dest="bash_file",help="run bash file to get status first",default="Off")
parser.add_option("-r","--resub",dest="resub_file",help="file where resub suggestions will be put",default="resubFile.sh")
(options,args) = parser.parse_args()

if options.bash_file!='Off':
  os.system("bash %s > %s"%(options.bash_file,options.input_file)) 

readFile = open(options.input_file)
resubFile = open(options.resub_file,'w')
resubFile.write("#!/bin/bash\n")

job_status={'Done':0,'Cleared':0,'Aborted':0,'Submitting':0,'Submitted':0,'Ready':0,'Running':0,'Pending':0}
job_action={'Aborted':0,'Cleared':0,'Terminated':0,'SubRequested':0,'SubSuccess':0,'SubFailed':0,'Cleaned':0,'KillSuccess':0,}
failedJobs=[]
njobs=0

pp=pprint.PrettyPrinter(indent=4)

for line in readFile.readlines():
  if "working directory" in line:
    longpath=(line.split()[2]).split('/')
    fileName=longpath[len(longpath)-2]
    print '-------------------------------'
    print fileName
    failedJobs=[]
    # reset counters
    njobs=0
    for key in job_status.keys():
      job_status[key]=0
  elif "ExitCodes Summary" in line:
    print 'Total Jobs: %4d'%njobs
    for key in sorted(job_status.iterkeys()):
      if job_status[key]>0:
        print '%6d jobs %s'%(job_status[key],key)
    print 'Failed Jobs: ', failedJobs
    recommend_action="crab -c %s/%s"%(fileName,fileName)
    if job_status['Submitting']>0 or job_status['Submitted']>0:
      recommend_action+=' -kill'
      for j in failedJobs: 
        recommend_action+=' %d,'%j
      recommend_action=recommend_action[:-1]
      resubFile.write(recommend_action+'\n') 
    recommend_action="crab -c %s/%s"%(fileName,fileName)
    if job_status['Aborted'] or job_status['Submitting']>0 or job_status['Submitted']>0:
      recommend_action+=' -resubmit'
      for j in failedJobs: recommend_action+=' %d,'%j
      recommend_action=recommend_action[:-1]
      resubFile.write(recommend_action+'\n') 
  else:
    jobInfo=line.split()
    if len(jobInfo)<1:
      continue
    if not jobInfo[0].isdigit():
      continue
    jobNum=jobInfo[0]
    njobs+=1
    if len(jobInfo)<6:
      exit_code=10000
    else:
      exit_code=jobInfo[5]
    if exit_code!="0" and exit_code!=" ":
      failedJobs.append(int(jobNum))
    for key in job_status.keys():
      if jobInfo[2]==key:
        job_status[key]+=1

print 'Recommended action written into resubFile.txt'
