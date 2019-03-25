#include <string>
#include <iostream>
#include <fstream>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooExtendPdf.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TChain.h"

using namespace RooFit ;
using namespace std ;

class HiggsMuMuFit{
public :

  HiggsMuMuFit(TString inputFileList="foo.txt", TString outFileName="histo.root", TString isbkgstr="F", vector<TString> selection={"cat_index==7"}, TString proc="ggH", int category=1);
  virtual ~HiggsMuMuFit();
  void     sigfit();
  void     sigfit_DCB();
  void     bkgfit();
  void     bkgfit_BWZReduxMBern();
  void     getpdf();
  void     set(TString inputFileList, TString outFileName, TString isbkgstr, vector<TString> selection, TString proc, int category);
  TChain* loader(const string& inFile_name, const string& chain_name);
private:
  string infile;
  TString outfile;
  vector<TString> cuts;
  TString isbkg;
  TString procname;
  int ncat;
};

HiggsMuMuFit::HiggsMuMuFit(TString inputFileList, TString outFileName, TString isbkgstr, vector<TString> selection, TString proc, int category)
{
  infile = inputFileList.Data();
  outfile = outFileName;
  cuts = selection;
  isbkg = isbkgstr;
  procname = proc;
  ncat = category; 
}

void HiggsMuMuFit::set(TString inputFileList, TString outFileName, TString isbkgstr, vector<TString> selection, TString proc, int category)
{
  infile = inputFileList;
  outfile = outFileName;
  cuts = selection;
  isbkg = isbkgstr;
  procname = proc;
  ncat = category;
}

TChain* HiggsMuMuFit::loader(const string& inFile_name, const string& chain_name)
{
    TChain* chain = new TChain(chain_name.c_str());
    TString in_name(inFile_name);
    if(in_name.Contains("root")) {
        chain->Add(inFile_name.c_str());
        cout << "total events: " << chain->GetEntries() << " in " << inFile_name.c_str() << endl;
        return chain;
    }
    fstream input(inFile_name.c_str(), fstream::in);
    string file_name;
    int ncounter = 0;
    while (input >> file_name){
        // cout << "adding: " << file_name << endl;
        chain->Add(file_name.c_str());
        ncounter ++;
    }
    cout << "total events: " << chain->GetEntries() << " in " << ncounter << " files." << endl;
    input.close();
    return chain;
}

HiggsMuMuFit::~HiggsMuMuFit()
{
}
