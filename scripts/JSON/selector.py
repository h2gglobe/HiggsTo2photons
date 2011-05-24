#! /usr/bin/env python

import sys, csv, json, time

'''
This function will take a json file with objects formated as 
[run1, {[lumi1_1_min,lumi1_1_max],...[lumi1_N_min,lumi1_N_max]}] 
and write out to funcText.txt a function which returns true if
the run and lumis are in the ranges given by the json file.
It is intended to be placed into an analysis cc.h.


'''

if(len(sys.argv)==1):
  print "Please enter the json file as an arg."
  sys.exit()

# json file input
jf = open(sys.argv[1], 'r')
my_jf_dict = json.load(jf)

runs = my_jf_dict.keys()
runs.sort()

output = 'allGoodLumi_'+str(time.strftime("%Y-%m-%d_%H%M%S"))+'.txt'
f = open(output,'w')
f.write('bool LoopAll::lumiRunSelection(Int_t jentry) {\n\n')
f.write('  b_run->GetEntry(jentry);\n')
f.write('  b_lumis->GetEntry(jentry);\n\n')

# loop over runs
#for run, lumis in my_jf_dict.items():
for run in runs:
  lumis = my_jf_dict[run]
  # loop over lumi section ranges
  for lumiRange in lumis:
    if type(lumiRange) == type([]) and len(lumiRange) == 2:
      f.write('  if (run == '+str(run)+' && lumis>='+str(lumiRange[0])+' && lumis<='+str(lumiRange[1])+') return true;\n')
    else:
      print "Something wrong with format at "+str(run)+":"+str(lumiRange)
      print "\ttype(lumiRange) is "+str(type(lumiRange))
      print "\tlen(lumiRange) is "+str(len(lumiRange))
      print "\tShould be type = list, len = 2"

f.write('  return false;\n}\n')
f.close()
print output+" is complete"

