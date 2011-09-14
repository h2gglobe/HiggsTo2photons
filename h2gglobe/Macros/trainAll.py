#!/usr/bin/env python

mass_list = [120] 

#ada paramameters optimising
ada_list = [.2,.4,.6,.8,1.0] 
tree_list =[128,256,384,512]
depth_list= [2,3,4,5,6]
prune_list = ["NoPruning","CostComplexity"]

from subprocess import call
import os
pwd = os.getcwd()

for mass in mass_list : 
    for tree in tree_list : 
        for depth in depth_list :
            for prune in prune_list : 
                for ada in ada_list : 
                    call(["qsub", "-q", "hepmedium.q", "subTMVAToBatch.sh", str(pwd) , str(mass), str(tree), str(depth), str(prune), str(ada)])

