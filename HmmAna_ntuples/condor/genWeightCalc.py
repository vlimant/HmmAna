import os,sys,re,fileinput,string,shutil
import numpy as np
import ROOT as rt

datasets = [#"ttTo2l2v_2018",
    #        "ttsl_2018",
   #         "ggH_2018",
  #          "ttH_2018",
 #           "WplusH_2018",
#            "WminusH_2018",
#"ttTo2l2v_2017",
 #           "ttsl_2017",
            "ggH_2017",
   #         "ttH_2017",
    #        "WplusH_2017",
     #       "WminusH_2017",
    #"ZH_2017",

#"ttTo2l2v_2016",
#            "ttsl_2016",
 #           "ggH_2016",
 #           "ttH_2016",
  #          "WplusH_2016",
  #          "WminusH_2016",
  #  "ZH_2016",

     #       "WZTo1L1Nu2Q",
#            "ZZ",
       #     "WZTo3LNu",
        #    "WZTo2L2Q",
         #   "WWTo2L2Nu",
            #"WWToLNuQQ",
          #  "WWW_4F",
           # "WWZ_4F",
  #         "ZH_2018",
 # "TTWJetsToLNu",
 # "TTZToLLNuNu"
 #   "DYJetsToLL_VBFfilter_2018",
 #   "DYJetsToLL_VBFfilter_2017",
 #   "DYJetsToLL_VBFfilter_2016",
#"VBFHToMuMu_2016",
#"VBFHToMuMu_2017",
#"VBFHToMuMu_2018"

]

for data in datasets:
    filepath = data+".txt"
    with open(filepath) as fp:  
        rootFiles = fp.readlines()
        sum=0
        for line in rootFiles:
            line = line.rstrip("\n")
            print line
            in_file = rt.TFile.Open(line)
            hist = in_file.Get("h_sumOfgenWeight")
            sum=sum+hist.Integral()
            in_file.Close()
        print data, sum
