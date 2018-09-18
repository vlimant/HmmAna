============================================

CONDOR scripts

===========================================

1. Before running condor, make sure to change the following items:

2 . give appropriate addresses in the run_myprog.sh, proto_condor_submit and makecondorsubmit.py

3.  When running over data, use the option "T" for IsData in makeCondorSubmit.py. While running over mc, use "F"

4. to initiate condor jobs, do python makeCondoSubmit.py

5. For every data[0] name in the makeCondorSubmit, make sure to have the appropriate runlist file in the condor directory and the latest version of analyzeHmm (executable)

6. Also create the following directories before submitting jobs

   a. condor/condor_output/condor_logs

   b. condor/condor_submit

7. Make sure to have the RoccoR2017.txt in condor/condor_output/condor_logs directory.

8. You can resubmit failed jobs by doing condor_submit condor_submit/submit_condor_job (appropriate job name)
