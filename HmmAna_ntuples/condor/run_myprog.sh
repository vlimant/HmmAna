#!/bin/bash

#set filelist=${1}
#set outfile=${2}
#set data=${3}

#I think you need to set your root environment. Best is to checkout a CMSSW package, cmsenv will give you access to root
#cd /home/spandey/t3store/public/CMSSW/CMSSW_7_5_0_pre5/src
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#eval `scram runtime -sh`
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc630
#ulimit -c 0
#eval `scram runtime -sh`
#echo `which root`
cd /data/idutta/CMSSW_9_4_9/src/HmmAna_ntuples/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630
ulimit -c 0
eval `scram runtime -sh`
echo `which root`

/data/idutta/CMSSW_9_4_9/src/HmmAna_ntuples/condor/analyzeHiggsMuMu $1 $2 $3 $4
