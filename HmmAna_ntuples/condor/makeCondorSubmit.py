import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
datasets = [["DYJetsToLL_VBFfilter_2017", "DYJetsToLL_VBFfilter_2017"]
            #["DYJetsToLL_VBFfilter_2018","DYJetsToLL_VBFfilter_2018"]
            #["DYJetsToLL_VBFfilter_2016","DYJetsToLL_VBFfilter_2016"]
            #["DYJetsToLL","DYJetsToLL"]
            #["ggH","ggH"]
            #["VBFHToMuMu","VBFHToMuMu"]
            #["ZZ","ZZ"]
            #["ttTosemileptonic","ttTosemileptonic"]
            #["ZH","ZH"]
            #["ttH","ttH"]
            #["WminusH","WminusH"]
            #["WplusH","WplusH"]
            #["ttTo2l2v","ttTo2l2v"]
            #["WZTo1L1Nu2Q","WZTo1L1Nu2Q"]
            #["WZTo3LNu","WZTo3LNu"]
            #["WZTo2L2Q","WZTo2L2Q"]
            #["WWTo2L2Nu","WWTo2L2Nu"]
            #["WWTo2L2Nu_Up","WWTo2L2Nu_Up"]
            #["WWToLNuQQ","WWToLNuQQ"]
            #["WWW_4F","WWW_4F"]
            #["WWZ_4F","WWZ_4F"]
            #["TTZToLLNuNu","TTZToLLNuNu"]
            #["TTWJetsToLNu","TTWJetsToLNu"]
]

NSections = 10
readFiles = ""

for data in datasets:
    jobidx = 0
    if ( data[0]=="DYJetsToLL_VBFfilter_2017"):
        dataname  = "DYJetsToLL_VBFfilter_2017"
        inputfname = "DYJetsToLL_VBFfilter_2017.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    
    elif ( data[0]=="DYJetsToLL_VBFfilter_2018"):
        dataname = "DYJetsToLL_VBFfilter_2018"
        inputfname = "DYJetsToLL_VBFfilter_2018.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1    

    elif ( data[0]=="DYJetsToLL_VBFfilter_2016"):
        dataname = "DYJetsToLL_VBFfilter_2016"
        inputfname = "DYJetsToLL_VBFfilter_2016.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 1

    elif ( data[0]=="DYJetsToLL"):
        dataname = "DYJetsToLL"
        inputfname = "DYJetsToLL.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ggH"):
        dataname = "Glu"
        inputfname = "ggH.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="VBFHToMuMu"):
        dataname = "VBFHToMuMu"
        inputfname = "VBFHToMuMu.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="WZ"):
        dataname = "WZ"
        inputfname = "WZ.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    elif ( data[0]=="ZZ"):
        dataname = "ZZ"
        inputfname = "ZZ.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
    
    elif ( data[0]=="ttTosemileptonic"):
        dataname = "ttTosemileptonic"
        inputfname = "ttTosemileptonic.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ZH"):
        dataname = "ZH"
        inputfname = "ZH.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ttH"):
        dataname = "ttH"
        inputfname = "ttH.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WminusH"):
        dataname = "WminusH"
        inputfname = "WminusH.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WplusH"):
        dataname = "WplusH"
        inputfname = "WplusH.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ttTo2l2v"):
        dataname = "ttTo2l2v"
        inputfname = "ttTo2l2v.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WZTo1L1Nu2Q"):
        dataname = "WZTo1L1Nu2Q"
        inputfname = "WZTo1L1Nu2Q.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WZTo3LNu"):
        dataname = "WZTo3LNu"
        inputfname = "WZTo3LNu.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="WZTo2L2Q"):
        dataname = "WZTo2L2Q"
        inputfname = "WZTo2L2Q.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="WWTo2L2Nu"):
        dataname = "WWTo2L2Nu"
        inputfname = "WWTo2L2Nu.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWTo2L2Nu_Up"):
        dataname = "WWTo2L2Nu_Up"
        inputfname = "WWTo2L2Nu_Up.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WWToLNuQQ"):
        dataname = "WWToLNuQQ"
        inputfname = "WWToLNuQQ.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWW_4F"):
        dataname = "WWW_4F"
        inputfname = "WWW_4F.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="WWZ_4F"):
        dataname = "WWZ_4F"
        inputfname = "WWZ_4F.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ttJets_DiLept"):
        dataname = "ttJets_DiLept"
        inputfname = "ttJets_DiLept.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="TTZToLLNuNu"):
        dataname = "TTZToLLNuNu"
        inputfname = "TTZToLLNuNu.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="TTWJetsToLNu"):
        dataname = "TTWJetsToLNu"
        inputfname = "TTWJetsToLNu.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10
  
    elif ( data[0]=="MuData1"):
        dataname = "MuData1"
        inputfname = "Data1_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="MuData2"):
        dataname = "MuData2"
        inputfname = "Data2_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="MuData3"):
        dataname = "MuData3"
        inputfname = "Data3_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData4"):
        dataname = "MuData4"
        inputfname = "Data4_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData5"):
        dataname = "MuData5"
        inputfname = "Data5_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData6"):
        dataname = "MuData6"
        inputfname = "Data6_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData7"):
        dataname = "MuData7"
        inputfname = "Data7_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="MuData8"):
        dataname = "MuData8"
        inputfname = "Data8_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="MuData9"):
        dataname = "MuData9"
        inputfname = "Data9_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="MuData10"):
        dataname = "MuData10"
        inputfname = "Data10_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData11"):
        dataname = "MuData11"
        inputfname = "Data11_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData12"):
        dataname = "MuData12"
        inputfname = "Data12_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData13"):
        dataname = "MuData13"
        inputfname = "Data13_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData14"):
        dataname = "MuData14"
        inputfname = "Data14_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="MuData15"):
        dataname = "MuData15"
        inputfname = "Data15_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="MuData16"):
        dataname = "MuData16"
        inputfname = "Data16_Muon.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal

    print "Dataset ",  data[0], " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    outDir="/mnt/hadoop/store/user/idutta/HmmAna/HistogramRootFiles/"+data[0]+"3May2019"
    print outDir
    #if not os.path.exists(outDir):
        #mkdirComm="scram unsetenv -sh;gfal-mkdir " +outDir
        #os.system(mkdirComm)
    outdir ="HmmAna/HistogramRootFiles/"+data[0]+"3May2019"
    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/data/idutta/CMSSW_9_4_9/src/HmmAna/master/HmmAna_ntuples/condor/condor_submit/runList_"+data[0]+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+data[0]+"_"+str(jobidx)
        jobName      = "3May2019"+data[0]+"_job"+str(jobidx)
        outHistFile = data[0]+"_job"+str(jobidx)+".root"
        isData       ="F"
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            line=line.replace("ROOTOUT",outHistFile)
            line=line.replace("DATANAME",dataname)
            line=line.replace("ISDATA",isData)
            line=line.replace("OUTDIR",outdir)
            print line.rstrip()
                        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
