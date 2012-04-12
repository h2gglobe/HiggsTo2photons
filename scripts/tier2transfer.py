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

myHadoop="/hadoop/cms/store/user/capalmer/h2g_V11_04_05_reloadred"

#transferList=["PhotonRun2011B_Clean_49.root"]


for dir in dirs:
    dir=dir.split("\n")[0]
    #if dir.find("_49") is -1:
    #    continue
    files = os.popen("/usr/bin/nsls " + sys.argv[1] + "/" + dir).readlines()
    filesinfo = os.popen("/usr/bin/nsls -l " + sys.argv[1] + "/" + dir).readlines()
    for fileinfo in filesinfo:
        fileinfo = fileinfo.split()
    

    MaxThreads=4*int(os.popen("cat /proc/cpuinfo | grep processor | awk '{print $3}' | tail -1 | xargs -i echo '{}+1' | bc").readlines()[0])
    sleepcount=0
    print MaxThreads, len(files)
    while len(filesinfo) > 0:
        threads = int(os.popen("ps | grep lcg-cp | wc | awk '{print $1}'").readlines()[0])
        if threads < MaxThreads:
            sleepcount=0
            fileinfo = filesinfo.pop().split()
            file = fileinfo[8]
            filesize = fileinfo[4]
            print file, threads, len(files)

            if file.find("_49") is -1:
                continue
            # check files existence and size
            line = """
lcg-ls -l -V cms -D srmv2 \"srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=%s/%s/%s\"
""" % (myHadoop, os.path.basename(sys.argv[1])+"/"+dir, file)
            input,checkOutput = os.popen2(line)
            hadoop_file_info = checkOutput.readlines()
            toCopy = False
            hadoop_file_size = 0 
            if (len(hadoop_file_info) == 0):
                toCopy = True
                print "Didn't find the file"
            else:
                hadoop_file_size = hadoop_file_info[0].split()[4]
                if hadoop_file_size != filesize:
                    print "hadoop_file_size, filesize  "+str(hadoop_file_size)+"  "+str(filesize)
                    toCopy = True
                    print "File is not the right size"
             
            
            if toCopy == False:
                continue
            if ("duplicate" in file):
                continue
            if ("empty" in file):
                continue
            line="""
lcg-cp -v -V cms -U srmv2 -T srmv2 -n 6 \"srm://srm-cms.cern.ch/srm/managerv2?SFN=%s/%s\" \"srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=%s/%s/%s\" &
""" % (sys.argv[1]+"/"+dir, file, myHadoop, os.path.basename(sys.argv[1])+"/"+dir, file)
            print line
            os.system(line)
        else:
            sleepcount=sleepcount+1
            #if sleepcount>20:
            #  print "flush lcg-cp's... I'm tired of sleeping - "+str(sleepcount)
            #  os.system("ps | grep lcg-cp | awk '{print 'kill -9 '$1}' | sh")
            #else:
            print "MAX THREADS - Sleeping for 30 sec - "+str(sleepcount)
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



