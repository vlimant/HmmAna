import os,sys,re,fileinput,string,shutil
import numpy as np
import ROOT as rt

datasets = ["ttTo2l2v",
#            "ttTosemileptonic",
            "ggH",
            "ttH",
            "WplusH",
            "WminusH",
            "WZTo1L1Nu2Q",
#            "ZZ",
            "WZTo3LNu",
            "WZTo2L2Q",
            "WWTo2L2Nu",
            "WWToLNuQQ",
            "WWW_4F",
            "WWZ_4F"
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
