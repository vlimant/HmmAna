This code train the MVA for the Hmm VH-leptonic category (WH(W->ev) process):

1) to train:
python VHLepTraining_el.py 

2)to add the discriminator output score to the ntuple:
python EvalBDT_VHLep_el.py /eos/cms/store/user/nlu/Hmm/categorization/Oct03_v1/merged/VBFH.root
