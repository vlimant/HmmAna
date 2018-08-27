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
