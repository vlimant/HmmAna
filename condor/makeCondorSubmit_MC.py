import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
#datanames = [
#["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"],
#["TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8"],
#["GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8"],
#["WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8"],
#["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"],
#["WZTo3LNu_0Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8"]
#]

#NSections = 10
#readFiles = ""

#for data in datanames:
#    jobidx = 0
#    dataname  = data[0]
#    inputfname = "list/" + data[0]+".txt"

datanames = ["nanoAODlist"]

for data in datanames:
 print data
 filepath = data+".txt"
 with open(filepath) as fp:
  rootFiles = fp.readlines()
  for dataname in rootFiles:
    print "dataname: ", dataname
    jobidx = 0
    dataname = dataname.rstrip("\n")

    inputfname = "nanoAODlist/" + dataname +".txt"
    with open(inputfname) as inputfile:
         readFiles = inputfile.readlines()
         print "len(readFiles)", len(readFiles)
    NSections = 1

    label = "11Mar19_v1"

    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal
    print "Dataset ", dataname, " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    outDir="/mnt/hadoop/store/user/nlu/Hmm/ntuple/2016/Nano14Dec2018/MC/Prod_v1/"+dataname+"_"+label
    print outDir
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/data/nlu/work/Hmm/CMSSW_9_4_9/src/HmmAna_2016/condor/condor_submit/runList_"+dataname+"_"+label+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+dataname+"_"+label+"_"+str(jobidx)
        jobName      = label+dataname+"_job"+str(jobidx)
        outHistFile = dataname+"__"+str(jobidx)+".root"
        #isData       ="T"
        isData       ="F"
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            line=line.replace("ROOTOUT",outHistFile)
            line=line.replace("DATANAME",dataname)
            line=line.replace("ISDATA",isData)
            line=line.replace("OUTDIR",outDir)
            print line.rstrip()
        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
