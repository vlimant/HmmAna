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
#include "NtupleVariables.h"
#include "TH1D.h"
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
  void     Categorization(const char *, const char *, float , float);
  void     GenInfo(const char *, const char *);
  void     BookHistogram(const char *);
  void     clearTreeVectors();
 
  TFile *oFile;
  TTree *cattree;

  UInt_t    run;
  UInt_t    lumi;
  ULong64_t event;
  float genWeight = -9999.;
  int cat_index = -9999.;
  float Higgs_mass = -9999.;
  float Higgs_pt = -9999.;
  float Higgs_eta = -9999.;
  float extralep_pfRelIso03 = -9999.;
  float extralep_pt = -9999.;
  float extralep_eta = -9999.;
  float dRlepHiggs = -9999.;
  float dEtamm = -9999.;
  float dRmm = -9999.;
  float dPhimm = -9999.;  
  float MET_phi = -9999.;
  float MET_pt = -9999.;
  std::vector<int>   *l1_index;
  std::vector<int>   *l2_index;

  //define histograms
  TH1D *catyield;
  TH1D *h_diMuon_mass_ggH;
  TH1D *h_diMuon_mass_VBF;
  TH1D *h_diMuon_mass_VHHad;
  TH1D *h_diMuon_mass_ZHll;
  TH1D *h_diMuon_mass_WHlv;
  TH1D *h_diMuon_mass_ttHLep;
  TH1D *h_diMuon_mass_ttHHad;  
  TH1D *h_dRlepH;
  TH1D *h_extralep1_pt;
  TH1D *h_extralep1_eta;
  
  TH1D *h_gen_diMuon_m;
  TH1D *h_gen_extralep;
  TH1D *h_gen_dRlepH;
  TH1D *h_gen_extralep1_pt;
  TH1D *h_gen_extralep1_eta;

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
  TH1D *h_diMuon_mass_SR;
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
  TH1D *h_Nbjet;
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
  cattree =  new TTree("cattree","cattree");
  l1_index = new std::vector<int>();
  l2_index= new std::vector<int>(); 
  cattree->Branch("run", &run,"run/i");
  cattree->Branch("lumi", &lumi, "lumi/i");
  cattree->Branch("event", &event,"event/l");
  cattree->Branch("cat_index", &cat_index, "cat_index/i");
  cattree->Branch("genWeight", &genWeight,"genWeight/F");
  cattree->Branch("Higgs_mass", &Higgs_mass,"Higgs_mass/F");
  cattree->Branch("Higgs_pt", &Higgs_pt,"Higgs_pt/F");
  cattree->Branch("Higgs_eta", &Higgs_eta,"Higgs_eta/F");
  cattree->Branch("extralep_pfRelIso03", &extralep_pfRelIso03, "extralep_pfRelIso03/F");
  cattree->Branch("extralep_pt", &extralep_pt,"extralep_pt/F");
  cattree->Branch("extralep_eta", &extralep_eta,"extralep_eta/F");
  cattree->Branch("dRlepHiggs", &dRlepHiggs,"dRlepHiggs/F");
  cattree->Branch("dRmm", &dRmm,"dRmm/F");
  cattree->Branch("dEtamm", &dEtamm,"dEtamm/F");
  cattree->Branch("dPhimm", &dPhimm,"dPhimm/F");  
  cattree->Branch("l1_index", "vector<int>", &l1_index);
  cattree->Branch("l2_index", "vector<int>", &l2_index);
  cattree->Branch("MET_phi", &MET_phi, "MET_phi/F");
  cattree->Branch("MET_pt", &MET_pt, "MET_pt/F");

  catyield = new TH1D("h_category_yield","h_category_yield",10,0,10);
  h_diMuon_mass_ggH = new TH1D("h_diMuon_mass_ggH","diMuon_mass_ggH",10,120,130);
  h_diMuon_mass_VBF = new TH1D("h_diMuon_mass_VBF","diMuon_mass_VBF",10,120,130);
  h_diMuon_mass_VHHad = new TH1D("h_diMuon_mass_VHHad","diMuon_mass_VHHad",10,120,130);
  h_diMuon_mass_ZHll = new TH1D("h_diMuon_mass_ZHll","diMuon_mass_ZHll",10,120,130);
  h_diMuon_mass_WHlv = new TH1D("h_diMuon_mass_WHlv","diMuon_mass_WHlv",10,120,130);
  h_diMuon_mass_ttHLep = new TH1D("h_diMuon_mass_ttHLep","diMuon_mass_ttHLep",10,120,130);
  h_diMuon_mass_ttHHad = new TH1D("h_diMuon_mass_ttHHad","diMuon_mass_ttHHad",10,120,130);
  h_dRlepH = new TH1D("h_dRlepH","deltaR(extra lepton and Higgs)",10,0,10);
  h_extralep1_pt = new TH1D("h_extralep1_pt","extralep1_pt",10,0,100);
  h_extralep1_eta = new TH1D("h_extralep1_eta","extralep1_eta",10,-5,5);

  h_gen_diMuon_m = new TH1D("h_gen_higgs_mass","Higgs mass",10,110,150);
  h_gen_extralep = new TH1D("h_gen_extralep","number of extra lepton",5,0,5);
  h_gen_dRlepH = new TH1D("h_gen_dRlepH","gen_deltaR(extra lepton and Higgs)",10,0,10);
  h_gen_extralep1_pt = new TH1D("h_gen_extralep1_pt","gen_extralep1_pt",10,0,100);
  h_gen_extralep1_eta = new TH1D("h_gen_extralep1_eta","gen_extralep1_eta",10,-5,5);

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
  h_diMuon_mass_SR=new TH1D("diMuon_mass_SR","Mass for dimuon system",70,60.,200.);
  h_diMuon_mass_110To120=new TH1D("diMuon_mass_110To120","Mass for dimuon system",5,110.,120.);
  h_diMuon_mass_130To150=new TH1D("diMuon_mass_130To150","Mass for dimuon system",10,130.,150.);
  h_diMuon_mass_110To150=new TH1D("diMuon_mass_110To150","Mass for dimuon system",20,110.,150.);

  h_j1pt=new TH1D("j1pt","P_{T} for leading jet",100,0.0,1000.);
  h_j1phi=new TH1D("j1phi","#phi for leading jet",32,-3.2,3.2);
  h_j1eta=new TH1D("j1eta","#eta for leading jet",40,-4.,4.);

  h_j2pt=new TH1D("j2pt","P_{T} for sub-leading jet",100,0.0,1000.);
  h_j2phi=new TH1D("j2phi","#phi for sub-leading jet",32,-3.2,3.2);
  h_j2eta=new TH1D("j2eta","#eta for sub-leading jet",40,-4.,4.);

  h_Njet=new TH1D("Njet","Number of jets",11,0.0,11.);
  h_Nbjet=new TH1D("Nbjet","Number of b-jets",11,0.0,11.);
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

void HiggsMuMu::clearTreeVectors(){
   run = -9999;
   lumi = -9999;
   event = -9999;
   genWeight = -9999.;
   cat_index = -9999.;
   Higgs_mass = -9999.;
   Higgs_pt = -9999.;
   Higgs_eta = -9999.;
   extralep_pfRelIso03 = -9999.;
   extralep_pt = -9999.;
   extralep_eta = -9999.;
   dRlepHiggs = -9999.;
   dEtamm = -9999.;
   dRmm = -9999.;
   dPhimm = -9999.;
   MET_phi = -9999.;
   MET_pt = -9999.;
   l1_index->clear();
   l2_index->clear();
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
