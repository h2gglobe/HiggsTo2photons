#!/usr/bin/env python
import sys

from subCatsCommon import *

if not options.minCats:
  print "No minCats defined: Using default (%d)"%(4)
  options.minCats = 4

if not options.maxCats:
  print "No maxCats defined: Using default (%d)"%(6)
  options.maxCats = 6

config = readDatfile(options.datfile)

print "Building directory structure and config files"
for j in range(options.minCats,options.maxCats+1):
  jobname = "%dCategories"%(j)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  print "\tJob: %s"%jobname

  os.system('mkdir -p %s'%jobdir)

  datfile='%s/%dCats_%s'%(jobdir,j,options.datfile)
  f = open('%s'%(datfile),'w')
  copyDatfile(config, f, ["wsfile", "nSpinCats", "catBoundaries"])
  f.write("wsfile=%s/Workspace.root\n"%(jobdir))
  f.write("nSpinCats=%d\n"%j)
  f.close()

  prepareJobArray(jobname, jobdir, datfile, options.njobs, options.toyspj)

print "Building workspaces"
for j in range(options.minCats,options.maxCats+1):
  jobname = "%dCategories"%(j)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  datfile='%s/%dCats_%s'%(jobdir,j,options.datfile)
  print "\tJob: %s"%jobname

  os.chdir(jobdir)
  processCommand = "%s/bin/diyBuildWorkspace %s &> buildWorkspace.log"%(startDir, datfile)
  if options.dryRun: print processCommand
  else: os.system(processCommand)
  os.chdir(startDir)

print "Submitting jobs"
for j in range(options.minCats,options.maxCats+1):
  jobname = "%dCategories"%(j)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  datfile='%s/%dCats_%s'%(jobdir,j,options.datfile)
  print "\tJob: %s"%jobname

  submitJobArray(jobname, jobdir, options.njobs, options.queue)

