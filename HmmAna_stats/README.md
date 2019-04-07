Preliminary script to study the bias in the VBF category using MC events.


1. To create inputs:
./run_VBFbias_test_v2.sh

2. setup Higgs combineTool

3. run datacard and extract bias

datacard/VBF_mass_sculpting_check/bias.sh

study bias on signal: fit to DY+VBF MC:
datacard/VBF_mass_sculpting_check/VBFcat0_bkg.txt

study bkg spurious signal
fit to DY MC:
datacard/VBF_mass_sculpting_check/VBFcat0.txt

For these two datacards, I am using the inputs from here: /eos/cms/store/user/nlu/Hmm/ML/output_test, which is created by step 1.
