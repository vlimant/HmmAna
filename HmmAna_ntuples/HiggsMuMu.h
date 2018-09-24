/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 10 04:03:58 2018 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: ../condor/condor_output/condor_logs/ElData16_job0.root
//////////////////////////////////////////////////////////

#ifndef HiggsMuMu_h
#define HiggsMuMu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
class HiggsMuMu:public NtupleVariables  {
public :

  HiggsMuMu(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *isData="F");
  virtual ~HiggsMuMu();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *);
  void     BookHistogram(const char *);
 
  TFile *oFile;
  //define histograms
  TH1D *h_mu1pt;
  TH1D *h_mu2pt;
  TH1D *h_mu1eta;
  TH1D *h_mu2eta;
  TH1D *h_mu1phi;
  TH1D *h_mu2phi;
  TH1D *h_mu1mu2dR;
  TH1D *h_mu1mu2dPhi;
  TH1D *h_diMuon_pt;
  TH1D *h_diMuon_eta;
  TH1D *h_diMuon_phi;
  TH1D *h_diMuon_mass;

  TH1D *h_diMuon_mass_110To120;
  TH1D *h_diMuon_mass_130To150;
  TH1D *h_diMuon_mass_110To150;
  TH1D *h_j1pt;
  TH1D *h_j1phi;
  TH1D *h_j1eta;

  TH1D *h_j2pt;
  TH1D *h_j2phi;
  TH1D *h_j2eta;

  TH1D *h_Njet;
  TH1D *h_j1j2dR;
  TH1D *h_j1j2dPhi;

  TH1D *h_dijet_pt;
  TH1D *h_dijet_eta;
  TH1D *h_dijet_phi;
  TH1D *h_Mjj;
  
  TH1D *h_MET_pt;
  TH1D *h_METphi;
  TH1D *h_MET_sumEt;
};

#endif

#ifdef HiggsMuMu_cxx

void HiggsMuMu::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;0

 

  oFile = new TFile(outFileName, "recreate");
  //oFile->mkdir("Cutflow");
  //oFile->cd("Cutflow");
  h_mu1pt=new TH1D("mu1pt","P_{T} for leading muon",100,0.0,1000.);
  h_mu1phi=new TH1D("mu1phi","#phi for leading muon",32,-3.2,3.2);
  h_mu1eta=new TH1D("mu1eta","#eta for leading muon",40,-4.,4.);

  h_mu2pt=new TH1D("mu2pt","P_{T} for sub-leading muon",100,0.0,1000.);
  h_mu2phi=new TH1D("mu2phi","#phi for sub-leading muon",32,-3.2,3.2);
  h_mu2eta=new TH1D("mu2eta","#eta for sub-leading muon",40,-4.,4.);
  h_mu1mu2dR=new TH1D("mu1mu2dR","#Delta R between two leading muons",15,0.,4.);
  h_mu1mu2dPhi=new TH1D("mu1mu2dPhi","#Delta #phi for two leading muons",32,-3.2,3.2);
  h_diMuon_pt=new TH1D("diMuon_pt","P_{T} for dimuon system",100,0.0,1000.);
  h_diMuon_eta=new TH1D("diMuon_phi","#phi for dimuon system",32,-3.2,3.2);
  h_diMuon_phi=new TH1D("diMuon_eta","#eta for dimuon system",40,-4.,4.);
  h_diMuon_mass=new TH1D("diMuon_mass","Mass for dimuon system",70,60.,200.);

  h_diMuon_mass_110To120=new TH1D("diMuon_mass_110To120","Mass for dimuon system",10,110.,120.);
  h_diMuon_mass_130To150=new TH1D("diMuon_mass_130To150","Mass for dimuon system",10,130.,150.);
  h_diMuon_mass_110To150=new TH1D("diMuon_mass_110To150","Mass for dimuon system",20,110.,150.);

  h_j1pt=new TH1D("j1pt","P_{T} for leading jet",100,0.0,1000.);
  h_j1phi=new TH1D("j1phi","#phi for leading jet",32,-3.2,3.2);
  h_j1eta=new TH1D("j1eta","#eta for leading jet",40,-4.,4.);

  h_j2pt=new TH1D("j2pt","P_{T} for sub-leading jet",100,0.0,1000.);
  h_j2phi=new TH1D("j2phi","#phi for sub-leading jet",32,-3.2,3.2);
  h_j2eta=new TH1D("j2eta","#eta for sub-leading jet",40,-4.,4.);

  h_Njet=new TH1D("Njet","Number of jets",11,0.0,11.);
  h_j1j2dR=new TH1D("j1j2dR","#Delta R between two leading jets",15,0.,4.);
  h_j1j2dPhi=new TH1D("j1j2dPhi","#Delta #phi for two leading jets",32,-3.2,3.2);

  h_dijet_pt=new TH1D("dijet_pt","P_{T} for dijet system",100,0.0,1000.);
  h_dijet_eta=new TH1D("dijet_phi","#phi for dijet system",32,-3.2,3.2);
  h_dijet_phi=new TH1D("dijet_eta","#eta for dijet system",40,-4.,4.);
  h_Mjj=new TH1D("diJet_mass","Mass for dijet system",100,60.,1000.);
  
  h_MET_pt=new TH1D("MET_pt","MET P_{T}",100,0.0,1000.);
  h_METphi=new TH1D("MET_phi","MET #phi",32,-3.2,3.2);
  h_MET_sumEt=new TH1D("MET_sumEt","MET Sum E_{T}",100,0.,1000.);
}

HiggsMuMu::HiggsMuMu(const TString &inputFileList, const char *outFileName, const char* dataset, const char *isData)
{
TChain *tree = new TChain("tree");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  NtupleVariables::Init(tree);
  BookHistogram(outFileName);
  }
Bool_t HiggsMuMu::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

HiggsMuMu::~HiggsMuMu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   oFile->Write();
   oFile->Close();
   
}

Long64_t HiggsMuMu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      // Notify();
   }
   return centry;
}
#endif // #ifdef HiggsMuMu_cxx
