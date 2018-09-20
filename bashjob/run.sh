#!/bin/bash
path='/afs/cern.ch/user/n/nlu/work/private/CMS/Hmm/CMSSW_9_4_9/src/HmmAna/master/bsub'
echo $path
cd $path
#echo "Starting job on " `date` #Date/time of start of job
#echo "Running on: `uname -a`" #Condor job is running on this node
#source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a tcsh script, use .csh instead of .sh
#export SCRAM_ARCH=slc6_amd64_gcc630
#eval 'scramv1 project CMSSW CMSSW_9_4_9' # This should not be needed. see [3]
#cd CMSSW_9_4_9/src/
#cmsenv
#voms-proxy-info -all
#voms-proxy-init -voms cms
#eval 'scramv1 runtime -sh' # cmsenv is an alias not on the workers # This also should not be needed, see [3]
#echo "CMSSW: "$CMSSW_BASE

#cd $path/workarea
input=$1
data=$2
path=$3
prefix=$4
opt=$5
echo $input $data $opt $path
#./sig_only_fit_data_perDetID_2018_test /eos/cms/store/group/dpg_hcal/comm_hcal/nlu/ntuples/Aug162018/MuonEGammaTOTEM_ECAL_nvtx_sel${ivtx}.root /eos/cms/store/group/dpg_hcal/comm_hcal/nlu/plots/tree_HcalCalHBHEMuonFilter-Aug16_MuonEGammaTOTEM/ outtree_${ieta}_${iphi}_${ivtx}.root $ieta $iphi $ivtx $ptcut > /eos/cms/store/group/dpg_hcal/comm_hcal/nlu/plots/tree_HcalCalHBHEMuonFilter-Aug16_MuonEGammaTOTEM/log_${ieta}_${iphi}_${ivtx}.log

#./analyzeHmm runList_data_muon_type.txt out_data_muon_type.root data T
#./analyzeHmm runList_MC.txt out_MC.root mc F
cp $input $path/
./analyzeHmm $input $path/${prefix}_out.root $data $opt > $path/${prefix}_log.txt
