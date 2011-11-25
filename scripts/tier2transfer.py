#! /usr/bin/python
import time
import os, sys

#print sys.argv[1]
#for recursion 
#dirs = os.popen("/usr/bin/nsls "+ sys.argv[1]+" | grep _newHLT").readlines()
dirs = os.popen("/usr/bin/nsls "+ sys.argv[1]).readlines()
print dirs
counter=1
NUMPARA=60
flag=0


for dir in dirs:
    dir=dir.split("\n")[0]
    if dir.find("WZH135") is -1:
      continue
    #if dir.find("2011B_newHLT") is not -1:
    #    continue
    #print dir,dir.find("WZH135")
    #if dir.find("WZH135") is 0 and flag is 0:
    #    flag=1
    #if flag is 1:
    files = os.popen("/usr/bin/nsls " + sys.argv[1] + "/" + dir).readlines()
    MaxThreads=3*int(os.popen("cat /proc/cpuinfo | grep processor | awk '{print $3}' | tail -1 | xargs -i echo '{}+1' | bc").readlines()[0])
    print MaxThreads, len(files)
    while len(files) > 0:
        threads = int(os.popen("ps | grep lcg-cp | wc | awk '{print $1}'").readlines()[0])
        if threads < MaxThreads:
            file = files.pop()
            file=file.split("\n")[0]
            print file, threads, len(files)
            if ("duplicate" in file):
                continue
            if ("empty" in file):
                continue
            line="""
lcg-cp -v -V cms -U srmv2 -T srmv2 -n 6 \"srm://srm-cms.cern.ch/srm/managerv2?SFN=%s/%s\" \"srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=/hadoop/cms/store/user/capalmer/h2g_V11_04_01/h2gred_nov18/%s/%s\" &
""" % (sys.argv[1]+"/"+dir, file, os.path.basename(sys.argv[1])+"/"+dir, file)
            os.system(line)
        else:
            print "MAX THREADS - Sleeping for 30 sec"
            time.sleep(30)

    #print files
    #for file in files:
    #    file=file.split("\n")[0]
    #    print dir,file
#   # for l in b.readlines():
    #    if ("duplicate" in file):
    #        continue
    #    if ("empty" in file):
    #        continue
    #    
    #    counter=counter+1
    #    line=""
    #    if counter%NUMPARA==0:
    #      line="""
#lcg-cp -v -V cms -U srmv2 -T srmv2 -n 6 \"srm://srm-cms.cern.ch/srm/managerv2?SFN=%s/%s\" \"srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=/hadoop/cms/store/user/capalmer/h2g_V11_04_01/%s/%s\" 
#""" % (sys.argv[1]+"/"+dir, file, os.path.basename(sys.argv[1])+"/"+dir, file)
#    #    else:
#    #      line="""
#lcg-cp -v -V cms -U srmv2 -T srmv2 -n 6 \"srm://srm-cms.cern.ch/srm/managerv2?SFN=%s/%s\" \"srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=/hadoop/cms/store/user/capalmer/h2g_V11_04_01/%s/%s\" &
#""" % (sys.argv[1]+"/"+dir, file, os.path.basename(sys.argv[1])+"/"+dir, file)
    #    os.system(line)
    #    if counter%NUMPARA==0 and counter>30:
    #        time.sleep(15)



