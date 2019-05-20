import os, sys
import ROOT as rt
from root_numpy import root2array, tree2array
#from root_pandas import read_root
import h5py 

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold

from xgboost import XGBClassifier
from xgboost import plot_tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_curve
import keras
import numpy as np
import numpy.lib.recfunctions as nlr
import pandas as pd

from matplotlib import pyplot as plt
import math
from math import pow

from keras.models import Sequential, Model
from keras.layers import Dense, Activation,Dropout, Input, BatchNormalization
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
try:
    import setGPU
except:
    os.environ['KERAS_BACKEND'] = 'tensorflow'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(sys.argv[3])


import advTrain_2018_defFunc

from keras.models import load_model
def output_file_nonTrain(y_pred_cls,y_pred_adv_cls,evt_weight, x_test_index,z_test,INPUT_FILE,proc_tree,alpha_str):
        # Create a lookup table for discriminator value by event number
        disc_lookup_signal = {}
        disc_lookup_signal_adv = {}
        proc_evt_weight={}
        sum=0
        h_sig_simple=rt.TH1F("simpleNN_VBF","simpleNN_VBF",100,0.,1.);
        h_sig_adv=rt.TH1F("advNN_VBF","advNN_VBF",100,0.,1.);
        for disc_val, disc_val_adv, evt_wt,mass, run, lumi, event in zip(y_pred_cls, y_pred_adv_cls, evt_weight,z_test, x_test_index[:,0],x_test_index[:,1], x_test_index[:,2]):
        
            disc_lookup_signal[(run, lumi, event)] = disc_val
            disc_lookup_signal_adv[(run, lumi, event)] = disc_val_adv
            proc_evt_weight[(run, lumi, event)] = evt_wt
            if (mass<120. or mass>130.):
                h_sig_simple.Fill(disc_val,evt_wt)
                h_sig_adv.Fill(disc_val_adv, evt_wt)
            sum+=evt_wt
        print "Sum of evt_Weight: ",sum
            # We write out the signal and background events to a new ROOT file.
            # For events in the test set, we append the GBM discriminator.
            # For events in the train set, we set the discriminator to -1
            # to avoid accidentally reusing the training set in the analysis.
        print len(disc_lookup_signal), len(disc_lookup_signal_adv)
        #print INPUT_FILE
        #print os.path.basename(INPUT_FILE)
        print alpha_str
        out_name = os.path.basename(INPUT_FILE).replace('.root', alpha_str+'_NNscore.root')
        out_proc = rt.TFile(out_name, 'RECREATE')
        out_proc_tree = proc_tree.CloneTree(0)
        
        
        # This is the boilerplate code for writing something
        # to a ROOT tree from Python.
        
        rt.gROOT.ProcessLine("struct MyStruct2{float disc_advNN;};")
        rt.gROOT.ProcessLine("struct MyStruct3{float disc_simpleNN;};")
        rt.gROOT.ProcessLine("struct MyStruct{double evt_wt;};")
        from ROOT import MyStruct3, MyStruct2, MyStruct
        s = MyStruct3()
        s_adv = MyStruct2()
        s_wt = MyStruct()
        disc_branch_proc = out_proc_tree.Branch('disc_simpleNN', rt.AddressOf(s, 'disc_simpleNN'), 'disc_simpleNN/F');
        disc_branch_sig_proc = out_proc_tree.Branch('disc_advNN', rt.AddressOf(s_adv, 'disc_advNN'), 'disc_advNN/F');
        evt_weight_branch_proc = out_proc_tree.Branch('evt_weight', rt.AddressOf(s_wt, 'evt_wt'), 'evt_weight/D');
            
        print "Writing new ROOT process file with discriminator appended"
        advTrain_2018_defFunc.fill_discriminator(proc_tree, out_proc_tree, disc_lookup_signal,disc_lookup_signal_adv, proc_evt_weight,s, s_adv,s_wt)
            
        # Cristian's code uses GetCurrentFile() for this part.
        # I will do that too just in case (paranoia).
        out_proc.cd()
        h_sig_simple.Write()
        h_sig_adv.Write()
        out_proc_tree.GetCurrentFile().Write()
        #signalNevents.Write()
        out_proc_tree.GetCurrentFile().Close()
            
        print 'done writing output'

def data_writer(fileName, x_mean,x_std):
    
    print fileName
    Sample_file = rt.TFile.Open(fileName)
    sample_tree = Sample_file.Get("cattree")
    datanp1 = advTrain_2018_defFunc.convert_data(sample_tree)
    data_frame1 = pd.DataFrame.from_records(datanp1)
    
    data_frame1.loc[:,variables] = (data_frame1.loc[:,variables]-x_mean)/x_std
    data1 = data_frame1.values
    print len(data1)
    data_reduced1 = data1[:,var_indices]
    data_mass1=data1[:,mass_indices]
    
    data_pred_simple = m.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv = m_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv1 = m1_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    data_pred_adv2 = m2_cls.predict(data_reduced1, batch_size=10000)[:, 0]
    #output_file_nonTrain(data_pred_simple,data_pred_adv,data[:,wt_var_indices], data[:,id_var_indices],Sample_file,sample_tree,'0p51em4')
    output_file_nonTrain(data_pred_simple,data_pred_adv1,data1[:,wt_var_indices], data1[:,id_var_indices],data_mass1,fileName,sample_tree,'1em5')
    output_file_nonTrain(data_pred_simple,data_pred_adv2,data1[:,wt_var_indices], data1[:,id_var_indices],data_mass1,fileName,sample_tree,'5em5')
    output_file_nonTrain(data_pred_simple,data_pred_adv,data1[:,wt_var_indices], data1[:,id_var_indices],data_mass1,fileName,sample_tree,'1em4')


m = load_model('simple_cls_DYnVBF.h5')
m_cls = load_model('advCls_1em4_DYnVBF.h5')
m1_cls = load_model('advCls_1em5_DYnVBF.h5')
m2_cls = load_model('advCls_5em5_DYnVBF.h5')

print(keras.__version__)
variables = ['softJet5','dRmm','dEtamm','dPhimm','M_jj','pt_jj','eta_jj','phi_jj','M_mmjj','pt_mmjj','eta_mmjj','phi_mmjj','dEta_jj','Zep','dRmin_mj', 'dRmax_mj'
                                   ,'dRmin_mmj','dRmax_mmj','dPhijj','leadingJet_pt','subleadingJet_pt',
                                   'leadingJet_eta','subleadingJet_eta','leadingJet_qgl','subleadingJet_qgl','cthetaCS','Higgs_pt','Higgs_eta']#,'Higgs_mass' ]
mass_var=['Higgs_mass']
id_variables = ['run','lumi','event']
wt_variables = ['genWeight']

SIGNAL_FILE='/bigdata/shared/idutta/VBFHToMuMu_2018.root'
BKG_FILE_2018='/bigdata/shared/idutta/DYJetsToLL_VBFfilter_2018.root'
sigfile = rt.TFile.Open(SIGNAL_FILE)
bkgfile_2018 = rt.TFile.Open(BKG_FILE_2018)
sigtree = sigfile.Get("cattree")
bkgtree3 = bkgfile_2018.Get("cattree")
signp =advTrain_2018_defFunc. convert(sigtree, 1)
bkgnp3 = advTrain_2018_defFunc.convert(bkgtree3, 0)

sig_frame0 = pd.DataFrame.from_records(signp)
bkg_frame3= pd.DataFrame.from_records(bkgnp3)

  
var_indices = [sig_frame0.columns.get_loc(v) for v in variables] # get positions of all the variables set above
mass_indices = [sig_frame0.columns.get_loc(v) for v in mass_var]
id_var_indices = [sig_frame0.columns.get_loc(v) for v in id_variables]
wt_var_indices = [sig_frame0.columns.get_loc(v) for v in wt_variables]

# Standardize only the actual variables not event ID or genweights
x_mean_sig0 = (sig_frame0.loc[:,variables]).mean()
x_std_sig0 = (sig_frame0.loc[:,variables]).std()

x_mean_bkg3 = (bkg_frame3.loc[:,variables]).mean()
x_std_bkg3 = (bkg_frame3.loc[:,variables]).std()


'''
inputfname='/nfshome/idutta/NN_sig_files_2.txt'
with open(inputfname) as inputfile:
    for line in inputfile:
        line=line.rstrip()
        print "Reading file:", line
        data_writer(line,x_mean_sig0,x_std_sig0)

inputfname='/nfshome/idutta/NN_bkg_files.txt'
with open(inputfname) as inputfile:
    for line in inputfile:
        line=line.rstrip()
        print line.rstrip()
        data_writer(line,x_mean_bkg3,x_std_bkg3)
'''
data_writer('/bigdata/shared/idutta/Data_2017.root',x_mean_bkg3,x_std_bkg3)
data_writer('/bigdata/shared/idutta/Data_2018.root',x_mean_bkg3,x_std_bkg3)
