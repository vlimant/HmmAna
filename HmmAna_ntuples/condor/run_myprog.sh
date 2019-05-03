#!/bin/sh
env 

cd /data/idutta/CMSSW_9_4_9/src/HmmAna/master/HmmAna_ntuples/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
#cd $TMP
/data/idutta/CMSSW_9_4_9/src/HmmAna/master/HmmAna_ntuples/analyzeHiggsMuMu $1 $TMP/$2 $3 $4

eval `scram unsetenv -sh`
gfal-copy -p -f --checksum-mode=both file://$TMP/$2 gsiftp://transfer.ultralight.org//store/user/idutta/$5/$2
