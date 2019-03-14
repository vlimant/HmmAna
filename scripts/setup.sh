cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src/
cmsenv
git cms-addpkg cms-nanoAOD:master-102X
git checkout -b nanoAOD cms-nanoAOD/master-102X
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/irenedutta23/HmmAna.git MyAnalysis/HmmAna
git clone https://gitlab.cern.ch/uhh-cmssw/fsr-photon-recovery.git FSRPhotonRecovery
