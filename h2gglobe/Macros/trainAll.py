#!/usr/bin/env python

mass_list = [120] 
tree_list = [128,256,512]
depth_list= [4,8,16,32]
prune_list = ["NoPruning"]#,"CostComplexity","ExpectedError"]
node_list = [4,8,16,32,64,128,1024]

#mass_list = [120] 
#tree_list = [128]
#depth_list= [4]
#prune_list = ["NoPruning"]#,"CostComplexity","ExpectedError"]
#node_list = [4,8]

from subprocess import call
import os
pwd = os.getcwd()

for mass in mass_list : 
    for tree in tree_list : 
        for depth in depth_list :
            for prune in prune_list : 
                for node in node_list : 
                    call(["qsub", "-q", "hepshort.q", "subTMVAToBatch.sh", str(pwd) , str(mass), str(tree), str(depth), str(prune), str(node)])

