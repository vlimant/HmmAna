#!/bin/bash
set +e
set +x

MINBIAS_XS=69200
MINBIAS_XS_UP=72660 #5% up
MINBIAS_XS_DOWN=65740 #5% down

GOLDEN_2016=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
PILEUP_2016=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt

GOLDEN_2017=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
PILEUP_2017=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt

GOLDEN_2018=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Rereco/Cert_314472-325175_13TeV_EarlyReReco2018ABC_Collisions18_JSON.txt
PILEUP_2018=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt

#compute the data PU profile with a given minbias cross section
pileup_single() {
    pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode true --minBiasXsec $3 --maxPileupBin 100 --numPileupBins 100 $4
}

#compute the variations and copy to a single file
pileup_variations() {
    pileup_single $1 $2 $MINBIAS_XS "nominal.root"
    pileup_single $1 $2 $MINBIAS_XS_UP "up.root"
    pileup_single $1 $2 $MINBIAS_XS_DOWN "down.root"
    rootcp up.root:pileup nominal.root:pileup_plus
    rootcp down.root:pileup nominal.root:pileup_minus
    mv nominal.root $3
}

pileup_variations $GOLDEN_2016 $PILEUP_2016 "RunII_2016_data.root"
pileup_variations $GOLDEN_2017 $PILEUP_2017 "RunII_2017_data.root"
pileup_variations $GOLDEN_2018 $PILEUP_2018 "RunII_2018_data.root"
