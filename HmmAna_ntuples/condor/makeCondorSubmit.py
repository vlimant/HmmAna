import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
datasets = [#["ttbar", "ttbar"]
            #["DYJetsToLL_M50_PT-100To250","DYJetsToLL_M50_PT-100To250"]
            #["DYJetsToLL_M50_PT-250To400","DYJetsToLL_M50_PT-250To400"]
            #["DYJetsToLL_M50_PT-400To650","DYJetsToLL_M50_PT-400To650"]
            #["DYJetsToLL_M50_PT-650ToInf","DYJetsToLL_M50_PT-650ToInf"]
            #["WW","WW"]
            #["WZ","WZ"]
            #["ZZ","ZZ"]
            #["ElData1","ElData1"]    
            #["ElData2","ElData2"]
            #["ElData3","ElData3"]
            #["ElData4","ElData4"]
            #["ElData5","ElData5"]
            #["ElData6","ElData6"]
            #["ElData7","ElData7"]
            #["ElData8","ElData8"]    
            #["ElData9","ElData9"]
            #["ElData10","ElData10"]
            #["ElData11","ElData11"]
            #["ElData12","ElData12"]
            #["ElData13","ElData13"]
            #["ElData14","ElData14"]
            #["ElData15","ElData15"]
            ["ElData16","ElData16"]
            #["MuData1","MuData1"]
            #["MuData2","MuData2"]
            #["MuData3","MuData3"]
            #["MuData4","MuData4"]
            #["MuData5","MuData5"]
            #["MuData6","MuData6"]
            #["MuData7","MuData7"]
            #["MuData8","MuData8"]
            #["MuData9","MuData9"]
            #["MuData10","MuData10"]
            #["MuData11","MuData11"]
            #["MuData12","MuData12"]
            #["MuData13","MuData13"]
            #["MuData14","MuData14"]
            #["MuData15","MuData15"]
            #["MuData16","MuData16"]
]

NSections = 10
readFiles = ""

for data in datasets:
    jobidx = 0
    if ( data[0]=="ttbar"):
        dataname  = "ttbar"
        inputfname = "ttbar.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 5

    
    elif ( data[0]=="DYJetsToLL_M50_PT-100To250"):
        dataname = "DYJetsToLL_M50_PT-100To250"
        inputfname = "DYJetsToLL_M-50_PT100to250_noDYNLO.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10    

    elif ( data[0]=="DYJetsToLL_M50_PT-250To400"):
        dataname = "DYJetsToLL_M50_PT-250To400"
        inputfname = "DYJetsToLL_M-50_PT250to400_noDYNLO.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="DYJetsToLL_M50_PT-400To650"):
        dataname = "DYJetsToLL_M50_PT-400To650"
        inputfname = "DYJetsToLL_M-50_PT400to650_noDYNLO.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="DYJetsToLL_M50_PT-650ToInf"):
        dataname = "DYJetsToLL_M50_PT-650ToInf"
        inputfname = "DYJetsToLL_M-50_PT650toInf_noDYNLO.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="WW"):
        dataname = "WW"
        inputfname = "WW.txt"
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
    
    elif ( data[0]=="ElData1"):
        dataname = "ElData1"
        inputfname = "Data1_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ElData2"):
        dataname = "ElData2"
        inputfname = "Data2_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ElData3"):
        dataname = "ElData3"
        inputfname = "Data3_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData4"):
        dataname = "ElData4"
        inputfname = "Data4_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData5"):
        dataname = "ElData5"
        inputfname = "Data5_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData6"):
        dataname = "ElData6"
        inputfname = "Data6_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData7"):
        dataname = "ElData7"
        inputfname = "Data7_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ElData8"):
        dataname = "ElData8"
        inputfname = "Data8_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

        
    elif ( data[0]=="ElData9"):
        dataname = "ElData9"
        inputfname = "Data9_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10



    elif ( data[0]=="ElData10"):
        dataname = "ElData10"
        inputfname = "Data10_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData11"):
        dataname = "ElData11"
        inputfname = "Data11_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData12"):
        dataname = "ElData12"
        inputfname = "Data12_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData13"):
        dataname = "ElData13"
        inputfname = "Data13_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData14"):
        dataname = "ElData14"
        inputfname = "Data14_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10

    elif ( data[0]=="ElData15"):
        dataname = "ElData15"
        inputfname = "Data15_Ele.txt"
        with open(inputfname) as inputfile:
            readFiles = inputfile.readlines()
            print "len(readFiles)", len(readFiles)
        NSections = 10


    elif ( data[0]=="ElData16"):
        dataname = "ElData16"
        inputfname = "runList.txt"
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

    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/data/idutta/CMSSW_9_4_9/src/HmmAna_ntuples/condor/condor_submit/runList_"+data[0]+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+data[0]+"_"+str(jobidx)
        jobName      = "9Sep2018"+data[0]+"_job"+str(jobidx)
        outHistFile = data[0]+"_job"+str(jobidx)+".root"
        isData       ="T"
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            line=line.replace("ROOTOUT",outHistFile)
            line=line.replace("DATANAME",dataname)
            line=line.replace("ISDATA",isData)
            print line.rstrip()
                        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
