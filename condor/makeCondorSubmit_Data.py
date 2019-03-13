import os,sys,re,fileinput,string,shutil

##             Dataset        Name   
datasets = [["Run2016D-17Jul2018-v1_Y03042019_v2r"],
            ["Run2016H-17Jul2018-v1_Y03042019_v2r"],
            ["Run2016C-17Jul2018-v1_Y03042019_v2r"],
            ["Run2016B-17Jul2018_ver1-v1_Y03042019_v2r"],
            ["Run2016B-17Jul2018_ver2-v1_Y03042019_v2r"],
            ["Run2016F-17Jul2018-v1_Y03042019_v2r"],
            ["Run2016G-17Jul2018-v1_Y03042019_v2r"],
            ["Run2016E-17Jul2018-v1_Y03042019_v2r"]
]

NSections = 10
readFiles = ""

for data in datasets:
    jobidx = 0
    dataname  = data[0]
    inputfname = "nanoAODlist/" + data[0]+".txt"
    with open(inputfname) as inputfile:
         readFiles = inputfile.readlines()
         print "len(readFiles)", len(readFiles)
    NSections = 1

    label = "07Feb19_v1"

    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal
    print "Dataset ",  data[0], " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    outDir="/mnt/hadoop/store/user/nlu/Hmm/ntuple/2016/Nano14Dec2018/Data/Prod_v1/"+data[0]+"_"+label
    print outDir
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = "/data/nlu/work/Hmm/CMSSW_9_4_9/src/HmmAna_2016/condor/condor_submit/runList_"+data[0]+"_"+label+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)
            #inputRunList.write("root://cmsxrootd.fnal.gov/"+line)

        condorSubmit = "condor_submit/submitCondor_"+data[0]+"_"+label+"_"+str(jobidx)
        jobName      = label+data[0]+"_job"+str(jobidx)
        outHistFile = data[0]+"__"+str(jobidx)+".root"
        isData       ="T"
        #isData       ="F"
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
