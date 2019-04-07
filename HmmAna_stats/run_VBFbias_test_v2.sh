#crate dataset for bkg
./FitHiggsMuMu DY.txt output_test DY BWZRedux,BWZRedux F 0.8
#crate dataset for Signal plus bkg
./FitHiggsMuMu VBFHDY.txt output_test VBFHDY BWZRedux,BWZRedux F 0.8
#create signal shape and yield using 2018 VBF signal MC
./FitHiggsMuMu /eos/cms/store/user/nlu/Hmm/ML/ntuples/VBFHToMuMu_2018_NNscore.root output_test VBFH BWZRedux,BWZRedux T 0.8
#create DY bkg norm using 2018 MC
./FitHiggsMuMu DY2018.txt output_test bkg2018BWZRedux BWZRedux,BWZRedux T 0.8
#create DY bkg shape using 2016+2017+2018 MC
./FitHiggsMuMu DY.txt output_test bkgBWZRedux BWZRedux,BWZRedux T 0.8
