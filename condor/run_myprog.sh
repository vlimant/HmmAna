#!/bin/bash

#based on https://caltech.teamwork.com/#notebooks/119930
#test using:
#TMP=`pwd` ./condor/run_myprog.sh /data/jpata/hmumu/CMSSW_10_2_10/src/MyAnalysis/HmmAna/condor/test.txt test.root DYJetsToLL_1 F

#Print out all bash commands
set -x

#Abort bash script on any error
set -e

echo "Looking inside job scratch dir: $TMP"
ls -al $TMP

BASE_DIR=/data/jpata/hmumu/CMSSW_10_2_10/src/MyAnalysis/HmmAna/
INPUT_FILELIST=$1
OUTPUT_FILENAME=$2
DATANAME=$3
ISDATA=$4
OUTDIR=$5

cd $BASE_DIR
source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc630
ulimit -c 0
eval `scram runtime -sh`
env 

echo "output filename: $TMP/$OUTPUT_FILENAME"

#Run NanoAOD postprocessing step
FILENAMES=`cat $TMP/$INPUT_FILELIST | xargs`
$BASE_DIR/scripts/nano_postproc.py --data_period 2017 --outdir $TMP $FILENAMES
ls -al

#Run private ntuple
$BASE_DIR/analyzeHmm $TMP/$INPUT_FILELIST $TMP/$OUTPUT_FILENAME $DATANAME $ISDATA
ls -al

#copy output to final location
mv $TMP/$OUTPUT_FILENAME $OUTDIR/
