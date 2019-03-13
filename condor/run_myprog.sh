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
cd /data/nlu/work/Hmm/CMSSW_9_4_9/src/HmmAna_2016/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630
ulimit -c 0
eval `scram runtime -sh`
env 
echo "$TMPDIR/$2"
/data/nlu/work/Hmm/CMSSW_9_4_9/src/HmmAna_2016/condor/analyzeHmm $1 $TMPDIR/$2 $3 $4

mv $TMPDIR/$2 $5
