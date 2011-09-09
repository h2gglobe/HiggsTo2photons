#!/bin/bash
# run with $1 = $PWD, $2 = .dat, $3 = number of jobs, $4 = this job id (0->njobs-1) 
hostname
echo "Started At -> "
date
echo "------------------------------------------------"
echo "Running MvaAnalysis from  $1 with dat file $2"
echo "Printing dat file"
cat $1/$2
source /vols/cms/grid/setup.sh
export X509_USER_PROXY=/home/hep/mj509/.myproxy
cd /vols/cms02/mjarvis/CMSSW_4_2_6
eval `scramv1 runtime -sh`
cd $1
python fitter.py -i $2 -n $3 -j $4

echo "Finished At -> "
date
echo "------------------------------------------------"
