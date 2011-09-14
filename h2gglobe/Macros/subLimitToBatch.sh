#!/bin/bash
# run with $1 = $PWD, $2 = mass, $3 = trees, $4 = depth, $5 = pruning, $6 = # adaboost
hostname
echo "Started At -> "
date
echo "------------------------------------------------"
source /vols/cms/grid/setup.sh
export X509_USER_PROXY=/home/hep/mj509/.myproxy
cd /vols/cms02/mjarvis/CMSSW_4_2_6
eval `scramv1 runtime -sh`
cd $1
python limit.py -M $2 -T $3 -D $4 -P $5 -A $6
echo "Finished At -> "
date
echo "------------------------------------------------"
