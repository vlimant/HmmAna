text2workspace.py VBFcat1.txt -o workspace_VBF_cat1.root
combine -n Obs -M MultiDimFit -m 125 workspace_VBF_cat1.root --algo=singles --robustFit=1 --setParameterRanges r=-50,50 -S 0
