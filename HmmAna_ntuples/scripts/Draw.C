#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "vector"
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"

void Draw(){
 
  TFile rtFile_ggH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ggH.root");
  TFile rtFile_VBFH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/VBFH.root");
  TFile rtFile_ZH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ZH.root");
  TFile rtFile_WmH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WminusH.root");
  TFile rtFile_WpH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WplusH.root");
  TFile rtFile_ttH("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ttH.root");
  TFile rtFile_DY("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/DY.root");
  TFile rtFile_ttl("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ttTosemileptonic.root");
  TFile rtFile_tt2l("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ttTo2l2v.root");
  TFile rtFile_WZTo1L1Nu2Q("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WZTo1L1Nu2Q.root");
  TFile rtFile_ZZ("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/ZZ.root");
  TFile rtFile_WZTo3LNu("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WZTo3LNu.root");
  TFile rtFile_WZTo2L2Q("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WZTo2L2Q.root");
  TFile rtFile_WWTo2L2Nu("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WWTo2L2Nu.root");
  TFile rtFile_WWToLNuQQ("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WWToLNuQQ.root");
  TFile rtFile_WWW_4F("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WWW_4F.root");
  TFile rtFile_WWZ_4F("/eos/cms/store/user/nlu/Hmm/categorization/Sept27/merged/WWZ_4F.root");

  float lumi = 41.529*1000.;
  TH1F* h_ggH = (TH1F*)rtFile_ggH.Get("h_category_yield");
  h_ggH->Scale(0.009605*lumi/217554238.5);

  TH1F* h_VBFH = (TH1F*)rtFile_VBFH.Get("h_category_yield");
  h_VBFH->Scale(0.000823*lumi/4506449.599577);
  
  TH1F* h_ZH = (TH1F*)rtFile_ZH.Get("h_category_yield");
  h_ZH->Scale(0.000192*lumi/234623.306747);

  TH1F* h_WpH = (TH1F*)rtFile_WpH.Get("h_category_yield");
  h_WpH->Scale(0.000183*lumi/259992.317749);

  TH1F* h_WmH = (TH1F*)rtFile_WmH.Get("h_category_yield");
  h_WmH->Scale(0.000116*lumi/162196.811523);

  TH1F* h_ttH = (TH1F*)rtFile_ttH.Get("h_category_yield");
  h_ttH->Scale(0.000110*lumi/155014.531738);

  TH1F* h_DY = (TH1F*)rtFile_DY.Get("h_category_yield");
  h_DY->Scale(5765.4*lumi/(3241270753030.957031+489144902631.884766));

  TH1F* h_ttl = (TH1F*)rtFile_ttl.Get("h_category_yield");
  h_ttl->Scale(6.871e+02*lumi/10876536248.0);

  TH1F* h_tt2l = (TH1F*)rtFile_tt2l.Get("h_category_yield");
  h_tt2l->Scale(85.656*lumi/623402174.0);

  TH1F* h_WZTo1L1Nu2Q = (TH1F*)rtFile_WZTo1L1Nu2Q.Get("h_category_yield");
  h_WZTo1L1Nu2Q->Scale(1.161e+01*lumi/352741934.219);

  TH1F* h_WZTo3LNu = (TH1F*)rtFile_WZTo3LNu.Get("h_category_yield");
  h_WZTo3LNu->Scale(4.42965*lumi/93694769.25);

  TH1F* h_WZTo2L2Q = (TH1F*)rtFile_WZTo2L2Q.Get("h_category_yield");
  h_WZTo2L2Q->Scale(5.595*lumi/255256973.25);

  TH1F* h_ZZ = (TH1F*)rtFile_ZZ.Get("h_category_yield");
  h_ZZ->Scale(16.523*lumi/1949768.0);
 
  TH1F* h_WWTo2L2Nu = (TH1F*)rtFile_WWTo2L2Nu.Get("h_category_yield");
  h_WWTo2L2Nu->Scale(12.46*lumi/177178.179688);

  TH1F* h_WWToLNuQQ = (TH1F*)rtFile_WWToLNuQQ.Get("h_category_yield");
  h_WWToLNuQQ->Scale(4.599e+01*lumi/405648754.016);

  TH1F* h_WWW_4F = (TH1F*)rtFile_WWW_4F.Get("h_category_yield");
  h_WWW_4F->Scale(0.2086*lumi/50039.244873);

  TH1F* h_WWZ_4F = (TH1F*)rtFile_WWZ_4F.Get("h_category_yield");
  h_WWZ_4F->Scale(0.1651*lumi/41205.3044434);
  
  TH1F* h_VVV = (TH1F*)h_WWZ_4F->Clone();
  h_VVV->Add(h_WWW_4F);

  TH1F* h_VV = (TH1F*)h_ZZ->Clone();
  h_VV->Add(h_WWToLNuQQ);
  h_VV->Add(h_WWTo2L2Nu);
  h_VV->Add(h_ZZ);
  h_VV->Add(h_WZTo2L2Q);
  h_VV->Add(h_WZTo1L1Nu2Q);
  h_VV->Add(h_WZTo3LNu);

  int cat = 3;  
  TString cat_name[] = {"ttHLep","ttHHad","ttHothers","ZHll","WHlv","VBF","VHHad","ggH"};

  Float_t vals[] = {(float)h_ggH->GetBinContent(cat),(float)h_VBFH->GetBinContent(cat),(float)h_ZH->GetBinContent(cat), (float)h_WpH->GetBinContent(cat),(float)h_WmH->GetBinContent(cat), (float)h_ttH->GetBinContent(cat),(float)h_DY->GetBinContent(cat), (float)h_ttl->GetBinContent(cat), (float)h_tt2l->GetBinContent(cat), (float)h_VV->GetBinContent(cat), (float)h_VVV->GetBinContent(cat)};
   Int_t colors[] = {1,2,3,4,5,6,7,8,9,10,11};
   Int_t nvals = sizeof(vals)/sizeof(vals[0]);

   TString name[] = {"ggH", "VBFH","ZH","WpH","WmH","ttH","DY","tt-semileptonic","tt2l2v","VV","VVV"};

   for(int i=0; i<11; i++){
      cout <<"expected number of events: "<<name[i]<<": "<<vals[i]<<endl;
   }
   TCanvas *cpie = new TCanvas("cpie","TPie test",700,700);
   TPie *pie1 = new TPie("pie1",
      "Pie with offset and no colors",nvals,vals,colors);
   //char lable[]={'ggH','VBFH','ZH','W-H','W+H','ttH','DY','ttl','tt2l','VV','VVV'};
   //pie1->SetLabels(lable);
   for (Int_t i=0;i<11;++i){
        pie1->GetSlice(i)->SetTitle(name[i]);
   }
   pie1->Draw("rsc");
   cpie->SaveAs("ttH"+cat_name[cat-1]+".png");
}
