To set up this code

cmsrel CMSSW_9_4_9

cd CMSSW_9_4_9/src

mkdir HmmAna

cd HmmAna

git clone git@github.com:irenedutta23/HmmAna master

make


To run the code

./analyzeHmm runList.txt out.root mc F

OR 

./analyzeHmm runList.txt out.root data T


============================================

CONDOR scripts

===========================================

Before running condor, make sure to change the following items:
1. give appropriate addresses in the run_myprog.sh, proto_condor_submit and makecondorsubmit.py
2. to initiate condor jobs, do python makeCondoSubmit.py
3. For every data[0] name in the makeCondorSubmit, make sure to have the appropriate runlist file in the condor directory and the latest version of analyzeHmm (executable)
4. Also create the following directories before submitting jobs

  a. condor/condor_output/condor_logs
  
  b. condor/condor_submit
  
5. You can resubmit failed jobs by doing condor_submit condor_submit/submit_condor_*job* (appropriate job name)
