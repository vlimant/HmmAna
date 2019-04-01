Preliminary script to study the bias in the VBF category using MC events.


1. To run:
run_VBFbias_test.sh

./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DY_VBF.root output_test_vbf_sb 2016_sbSExp SExp,SExp T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DY_VBF.root output_test_vbf_sb 2016_sbBWZGamma BWZGamma,BWZGamma T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DY_VBF.root output_test_vbf_sb 2016_sbBWZRedux BWZRedux,BWZRedux T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DY_VBF.root output_test_vbf_sb 2016_sbBWZ BWZ,BWZ T

These are running s+b fit to VBF+DY MC sample (DY_VBF.root), using bkg function SExp,BWZGamma,BWZRedux and BWZ functions.


./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DYJetsToLL_VBFfilter_2018_NNscore.root output_test_vbf_b 2016_sbSExp SExp,SExp T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DYJetsToLL_VBFfilter_2018_NNscore.root output_test_vbf_b 2016_sbBWZGamma BWZGamma,BWZGamma T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DYJetsToLL_VBFfilter_2018_NNscore.root output_test_vbf_b 2016_sbBWZRedux BWZRedux,BWZRedux T
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/VBF/categorization/DYJetsToLL_VBFfilter_2018_NNscore.root output_test_vbf_b 2016_sbBWZ BWZ,BWZ T

These are running s+b fit to DY MC sample, using bkg function SExp,BWZGamma,BWZRedux and BWZ functions.

output root file in output_test_vbf_b and output_test_vbf_sb contains the signal stength mu and the uncertainty from the fit (minos error). For example: output_test_vbf_b/2016_sbBWZRedux_result.root tree: result
signal stength mu from the fit: mu
lower error on mu from minos: errhi
upper error on mu from minos: errlo
err = 0.5*(errh+errlo)
