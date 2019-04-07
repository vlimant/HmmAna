#text2workspace.py VBFcat0.txt -o workspace_VBF_cat0.root
#combine -n Obs -M MultiDimFit -m 125 workspace_VBF_cat0.root --algo=singles --robustFit=1 --setParameterRanges r=-1,1.5 -S 0

text2workspace.py VBFcat0_bkg.txt -o workspace_VBF_cat0_bkg.root
combine -n Obs -M MultiDimFit -m 125 workspace_VBF_cat0_bkg.root --algo=singles --robustFit=1 --setParameterRanges r=-1,1.5 -S 0
#combine -n Obs -M MultiDimFit -m 125 workspace_VBF_cat0.root --algo=singles --robustFit=1 --setRobustFitStrategy 0 --setRobustFitTolerance 1.0 --setParameterRanges r=-50,50 -S 0 -v 2
