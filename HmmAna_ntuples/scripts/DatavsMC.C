#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"
using std::cout;
using std::endl;
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) ;
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame );

void DatavsMC()
{

  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  //lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_13TeV  = "41.529 fb^{-1}";
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=11; //11 = default: left-aligned
  //int iPos=33; // right-aligned
   //int iPos=22; // center-aligned
  //int iPos=0;//out of frame (in exceptional cases)
  //First declare the file names

    // references for T, B, L, R
  /* float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
  */
  TString f[77];
  f[0] = "DYJetsToLL_VBFfilter_2017.root"; //Fill in the appropriate names.
  f[1] = "VBFHToMuMu_2017.root";
  f[2] = "Data_2017.root";
  f[3] = "ttsl_2017.root";
  f[4] = "ttTo2l2v_2017.root";
  f[5] = "ggH_2017.root";
  f[6] = "WplusH_2017.root";
  f[7] = "WminusH_2017.root";
  f[8] = "ZH_2017.root";
  f[9] = "ttH_2017.root";
  f[10] = "ZZTo4L_2017.root";
  f[11]="ZZTo2L2Nu_2017.root";
  f[12] = "WZTo3LNu_2017.root";
  f[13]="WZTo2L2Q_2017.root";
  f[14]="WWTo2L2Nu_2017.root";
  f[15]="WWToLNuQQ_2017.root";
  f[16] = "WWW_4F_2017.root";
  f[17] = "ZZTo2L2Q_2017.root";
  f[18] = "TTWJetsToLNu_2017.root";
  f[19] = "TTZToLLNuNu_2017.root";
  f[20] = "WZZ_2017.root";
  f[21] = "ZZZ_2017.root";
  f[22] = "ttJets_DiLept_2017.root";
  f[23] = "EWK_2017.root";
  /* f[22] = "TTPrime_1200_JECup_Ht500_tZ.root";
  f[23] = "DYJetsToLL_M50_PT-250To400_JECup_Ht500.root";
  f[24] = "DYJetsToLL_M50_PT-400To650_JECup_Ht500.root";
  f[25] = "DYJetsToLL_M50_PT-650ToInf_JECup_Ht500.root";
  f[26] = "TTPrime_1100_JECup_Ht500_tZ.root";
  f[27] = "TTPrime_1300_JECup_Ht500_tZ.root";
  f[28] = "TTPrime_1400_JECup_Ht500_tZ.root";
  f[29] = "WW_JECup_Ht500.root";
  f[30] = "WZ_JECup_Ht500.root";
  f[31] = "ZZ_JECup_Ht500.root";

  f[32] = "DYJetsToLL_M50_PT-100To250_JECdown_Ht500.root"; //Fill in the appropriate names.
  f[33] = "ttbar_JECdown_Ht500.root";
  f[34] = "TTPrime_800_JECdown_Ht500_tZ.root";
  f[35] = "TTPrime_900_JECdown_Ht500_tZ.root";
  f[36] = "TTPrime_1000_JECdown_Ht500_tZ.root";
  f[37] = "TTPrime_1200_JECdown_Ht500_tZ.root";
  f[38] = "DYJetsToLL_M50_PT-250To400_JECdown_Ht500.root";
  f[39] = "DYJetsToLL_M50_PT-400To650_JECdown_Ht500.root";
  f[40] = "DYJetsToLL_M50_PT-650ToInf_JECdown_Ht500.root";
  f[41] = "TTPrime_1100_JECdown_Ht500_tZ.root";
  f[42] = "TTPrime_1300_JECdown_Ht500_tZ.root";
  f[43] = "TTPrime_1400_JECdown_Ht500_tZ.root";
  f[44] = "WW_JECdown_Ht500.root";
  f[45] = "WZ_JECdown_Ht500.root";
  f[46] = "ZZ_JECdown_Ht500.root";

   f[47] = "DYJetsToLL_M50_PT-100To250_JERup_Ht500.root"; //Fill in the appropriate names.
  f[48] = "ttbar_JERup_Ht500.root";
  f[49] = "TTPrime_800_JERup_Ht500_tZ.root";
  f[50] = "TTPrime_900_JERup_Ht500_tZ.root";
  f[51] = "TTPrime_1000_JERup_Ht500_tZ.root";
  f[52] = "TTPrime_1200_JERup_Ht500_tZ.root";
  f[53] = "DYJetsToLL_M50_PT-250To400_JERup_Ht500.root";
  f[54] = "DYJetsToLL_M50_PT-400To650_JERup_Ht500.root";
  f[55] = "DYJetsToLL_M50_PT-650ToInf_JERup_Ht500.root";
  f[56] = "TTPrime_1100_JERup_Ht500_tZ.root";
  f[57] = "TTPrime_1300_JERup_Ht500_tZ.root";
  f[58] = "TTPrime_1400_JERup_Ht500_tZ.root";
  f[59] = "WW_JERup_Ht500.root";
  f[60] = "WZ_JERup_Ht500.root";
  f[61] = "ZZ_JERup_Ht500.root";

  f[62] = "DYJetsToLL_M50_PT-100To250_JERdown_Ht500.root"; //Fill in the appropriate names.
  f[63] = "ttbar_JERdown_Ht500.root";
  f[64] = "TTPrime_800_JERdown_Ht500_tZ.root";
  f[65] = "TTPrime_900_JERdown_Ht500_tZ.root";
  f[66] = "TTPrime_1000_JERdown_Ht500_tZ.root";
  f[67] = "TTPrime_1200_JERdown_Ht500_tZ.root";
  f[68] = "DYJetsToLL_M50_PT-250To400_JERdown_Ht500.root";
  f[69] = "DYJetsToLL_M50_PT-400To650_JERdown_Ht500.root";
  f[70] = "DYJetsToLL_M50_PT-650ToInf_JERdown_Ht500.root";
  f[71] = "TTPrime_1100_JERdown_Ht500_tZ.root";
  f[72] = "TTPrime_1300_JERdown_Ht500_tZ.root";
  f[73] = "TTPrime_1400_JERdown_Ht500_tZ.root";
  f[74] = "WW_JERdown_Ht500.root";
  f[75] = "WZ_JERdown_Ht500.root";
  f[76] = "ZZ_JERdown_Ht500.root";
  
  */
  //Declare other constants, strings that you might need here.

  //Now let us define the plot name to overlay
  TString plotname[33];
  plotname[0]= "dijet_pt_VBF";
  plotname[1] = "dijet_eta_VBF";
  plotname[2] = "dijet_phi_VBF";
  plotname[3] = "mu1mu2dEta_VBF";
  plotname[4] ="mu1mu2dR_VBF";
  plotname[5] = "mu1mu2dPhi_VBF";
  plotname[6]= "diMuon_pt_VBF";
  plotname[7] = "diMuon_eta_VBF";
  plotname[8] = "diMuon_phi_VBF";
  plotname[9] = "h_diMuon_mass_VBF";
  plotname[10] = "h_diMuon_mass_110To150_VBF";
  plotname[11] = "diJet_mass_VBF";
  plotname[12] = "dijet_dEta_VBF";
  plotname[13] = "dijet_dPhijj_VBF";
  plotname[14] = "jet1_pt_VBF";
  plotname[15] = "jet2_pt_VBF";
  plotname[16] = "jet1_eta_VBF";
  plotname[17] = "jet2_eta_VBF";
  plotname[18] = "cthetaCS_VBF";
  plotname[19] = "dRmin_mj_VBF";
  plotname[20] = "dRmax_mj_VBF";
  plotname[21] = "dRmin_mmj_VBF";
  plotname[22] = "dRmax_mmj_VBF";
  plotname[23] = "Zep_VBF";
  plotname[24] = "di_muon_jet_mass_VBF";
  plotname[25] = "di_muon_jet_pt_VBF";
  plotname[26] = "di_muon_jet_eta_VBF";
  plotname[27] = "di_muon_jet_phi_VBF";
  plotname[28] = "jet1_qgl_VBF";
  plotname[29] = "jet2_qgl_VBF";
  plotname[30] = "softJet5_VBF";
  /*plotname[32]="Cutflow_DvsMC";*/
  gStyle->SetOptStat(0);
  //Also give fancy name for the axis titles
  TString xtitle= "P_{T} (GeV)";
  TString xtitle1= "#eta";
  TString xtitle2= "#phi";
  TString xtitle3= "H_{T} (GeV)";
  TString ytitle = "Number of Events"; // Or "Events"
  
  //Now let us open the files
  TFile *file[77];
  for(int i=0;i<24;i++) file[i] = new TFile(f[i]);
  
  //change directory if required.
  TDirectory *directory[77];
  for(int i=0;i<24;i++) directory[i]= (TDirectory*) file[i]->Get("VBF");
  // cout<<directory[0]<<endl;

  Double_t xbins[35]={0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300,400,500,700,1000.};
  
  TString name_unc;
  double sum=0.0;
  char hist_sum[20];
  TLegend *lg[32];
  char name[100];
  THStack *hs[32];
  TCanvas *c[32];
  TH1F *vFrame[32];
  TPad *pad1[32];
  TPad *pad2[32];
  TPad *pad3[32];
  TH1F *hd[32];
  TLine *l[32];
  //Now open the respective histograms from the files
  TH1F *h1[32],*h2[32],*h3[32],*h4[32];
  TH1F *h5[32],*h6[32],*h7[32],*h8[32];
  TH1F *h9[32],*h10[32],*h11[32],*h12[32];
  TH1F *h13[32],*h14[32],*h15[32],*h16[32];
  TH1F *h17[32],*h18[32],*h19[32],*h20[32];
  TH1F *h21[32],*h22[32],*h23[32],*h24[32];
  TH1F *h[32], *h_err[32], *h_temp;
  /*
  TH1F *h1_Up[32],*h1_Down[32],*h2_Up[32],*h2_Down[32];
  TH1F *h8_Up[32],*h8_Down[32],*h9_Up[32],*h9_Down[32];
  TH1F *h10_Up[32],*h10_Down[32],*h14_Up[32],*h14_Down[32];
  TH1F *h15_Up[32],*h15_Down[32],*h16_Up[32],*h16_Down[32];


  TH1F *h1_JECUp[32],*h1_JECDown[32],*h2_JECUp[32],*h2_JECDown[32];
  TH1F *h8_JECUp[32],*h8_JECDown[32],*h9_JECUp[32],*h9_JECDown[32];
  TH1F *h10_JECUp[32],*h10_JECDown[32],*h14_JECUp[32],*h14_JECDown[32];
  TH1F *h15_JECUp[32],*h15_JECDown[32],*h16_JECUp[32],*h16_JECDown[32];

  TH1F *h1_JERUp[32],*h1_JERDown[32],*h2_JERUp[32],*h2_JERDown[32];
  TH1F *h8_JERUp[32],*h8_JERDown[32],*h9_JERUp[32],*h9_JERDown[32];
  TH1F *h10_JERUp[32],*h10_JERDown[32],*h14_JERUp[32],*h14_JERDown[32];
  TH1F *h15_JERUp[32],*h15_JERDown[32],*h16_JERUp[32],*h16_JERDown[32];
  */
  double bin_err, stat_err,jec_err,jer_err;
  
  for(int j=0;j<31;j++){
    h1[j] = (TH1F*)directory[0]->Get(plotname[j]);
    h1[j]->Scale(1.); 
    sprintf(name,"h%i",j);
    h[j] = (TH1F*)h1[j]->Clone(name); //For MC // doing it here so that I can adjust the color of h[j]
    decorate(h1[j],"",ytitle,"",kYellow,2,kYellow,20,1);
   
  }

  for(int j=0;j<31;j++){
    h2[j] = (TH1F*)directory[1]->Get(plotname[j]);
    decorate(h2[j],"",ytitle,"",kViolet+1,2,kViolet+1,21,0);
  }

  for(int j=0;j<31;j++){
    h3[j] = (TH1F*)directory[2]->Get(plotname[j]);
    decorate(h3[j],"",ytitle,"",kBlack,2,kBlack,8,0);
  }
  
  for(int j=0;j<31;j++){
   h4[j] = (TH1F*)directory[3]->Get(plotname[j]);
    decorate(h4[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);
   }
 
  for(int j=0;j<31;j++){
     
    h5[j] = (TH1F*)directory[4]->Get(plotname[j]);
    decorate(h5[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);
     }
  for(int j=0;j<31;j++){
     
    h23[j] = (TH1F*)directory[22]->Get(plotname[j]);
    decorate(h23[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);
   }
  
  for(int j=0;j<31;j++){
   h6[j] = (TH1F*)directory[5]->Get(plotname[j]);
    decorate(h6[j],"",ytitle,"",kCyan+4,2,kCyan+4,21,0);
  }

  for(int j=0;j<31;j++){
    h7[j] = (TH1F*)directory[6]->Get(plotname[j]);
     decorate(h7[j],"",ytitle,"",kBlue,2,kBlue,21,0);
  }

  for(int j=0;j<31;j++){
    h8[j] = (TH1F*)directory[7]->Get(plotname[j]);
    decorate(h8[j],"",ytitle,"",kPink-7,2,kPink-7,21,0);
  }
    
  
  for(int j=0;j<31;j++){
    h9[j] = (TH1F*)directory[8]->Get(plotname[j]);
    decorate(h9[j],"",ytitle,"",kRed-2,2,kRed-2,21,0);
   
  }
  
  for(int j=0;j<31;j++){
   h10[j] = (TH1F*)directory[9]->Get(plotname[j]);
    decorate(h10[j],"",ytitle,"",kRed,2,kRed,21,0);
  }

  
  for(int j=0;j<31;j++){
    h11[j] = (TH1F*)directory[10]->Get(plotname[j]);
    decorate(h11[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }
    for(int j=0;j<31;j++){
    h12[j] = (TH1F*)directory[11]->Get(plotname[j]);
    decorate(h12[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    }
  
  for(int j=0;j<31;j++){
    h13[j] = (TH1F*)directory[12]->Get(plotname[j]);
    decorate(h13[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }

for(int j=0;j<31;j++){
    h14[j] = (TH1F*)directory[13]->Get(plotname[j]);
    decorate(h14[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }
  for(int j=0;j<31;j++){
    h15[j] = (TH1F*)directory[14]->Get(plotname[j]);
    decorate(h15[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }
  
  for(int j=0;j<31;j++){
   h16[j] = (TH1F*)directory[15]->Get(plotname[j]);
    decorate(h16[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }
  for(int j=0;j<31;j++){
    h17[j] = (TH1F*)directory[16]->Get(plotname[j]);
    decorate(h17[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
  }

  for(int j=0;j<31;j++){
    h18[j] = (TH1F*)directory[17]->Get(plotname[j]);
    decorate(h18[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
  }
  for(int j=0;j<31;j++){
   h19[j] = (TH1F*)directory[18]->Get(plotname[j]);
   decorate(h19[j],"",ytitle,"",kOrange,2,kOrange,21,1);
  }

  for(int j=0;j<31;j++){
    h20[j] = (TH1F*)directory[19]->Get(plotname[j]);
    decorate(h20[j],"",ytitle,"",kOrange,2,kOrange,21,1);
  }
  
  for(int j=0;j<31;j++){
     h21[j] = (TH1F*)directory[20]->Get(plotname[j]);
    decorate(h21[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
  }
  for(int j=0;j<31;j++){
    h22[j] = (TH1F*)directory[21]->Get(plotname[j]);
    decorate(h22[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
  }
  for(int j=0;j<31;j++){
    h24[j] = (TH1F*)directory[23]->Get(plotname[j]);
    decorate(h24[j],"",ytitle,"",kCyan,2,kCyan,21,1);
  }
  TLegend *lg1 = new TLegend(0.35,0.75,0.85,0.9);
  lg1->SetTextFont(132);
  lg1->SetBorderSize(0);
  lg1->SetNColumns(3);
  gStyle->SetLegendTextSize(0.04);
  lg1->AddEntry(h3[3],"Data","p");
  
  lg1->AddEntry(h1[3],"DY+jets","f");
  lg1->AddEntry(h24[3],"EWK LL+jets","f");
  lg1->AddEntry(h23[3],"t#bar{t}","f");
  lg1->AddEntry(h11[3],"VV ","f");
  
  lg1->AddEntry(h17[3],"VVV","f");
  
  lg1->AddEntry(h19[3],"ttX+Jets","f");
  

  lg1->AddEntry(h2[3],"VBF","l");
  lg1->AddEntry(h6[3],"ggH","l");
  lg1->AddEntry(h7[3],"W^{+}H","l");
  lg1->AddEntry(h8[3],"W^{-}H","l");
  lg1->AddEntry(h10[3],"ttH","l");
  lg1->AddEntry(h9[3],"ZH","l");
  for(int j=0;j<31;j++){
   
    sprintf(name,"hs%i",j);
    hs[j] = new THStack(name,plotname[j]); // Stack of all backgrounds
   
    
   
   
    hs[j]->Add(h19[j]);//ttw
    hs[j]->Add(h20[j]); //ttz
    hs[j]->Add(h17[j]); //www
    hs[j]->Add(h21[j]); //wzz
    hs[j]->Add(h22[j]); //zzz
   
   
    hs[j]->Add(h16[j]); //wwtolnuqq
    hs[j]->Add(h11[j]); // zz to 4l
     hs[j]->Add(h12[j]); //zz to 2l2v
    hs[j]->Add(h13[j]);//wz to 3l nu
    hs[j]->Add(h18[j]);//zz to 2l2q
    hs[j]->Add(h14[j]); //wz to 2l2q
   
    hs[j]->Add(h15[j]); //ww to 2l2v
    hs[j]->Add(h24[j]);//EWK
    hs[j]->Add(h23[j]);//ttJets dilept
    hs[j]->Add(h4[j]); //ttsl
   
    hs[j]->Add(h5[j]); //ttbar
    
    
    
    hs[j]->Add(h1[j]);
    

    sprintf(name,"c%i",j);
    c[j] = new TCanvas(name,name,800,800);
    /* c[j]->SetFillColor(0);
    c[j]->SetBorderMode(0);
    c[j]->SetFrameFillStyle(0);
    c[j]->SetFrameBorderMode(0);
    c[j]->SetLeftMargin( L/W );
    c[j]->SetRightMargin( R/W );
    c[j]->SetTopMargin( T/H );
    c[j]->SetBottomMargin( B/H );
    c[j]->SetTickx(0);
    c[j]->SetTicky(0);*/
    sprintf(name,"pad1%i",j);
    pad1[j]= new TPad(name, name, 0, 0.3, 1.0, 1.0);
    pad1[j]->SetBottomMargin(0.05); // Upper and lower plot are joined
    //pad1[j]->SetGridx();         // Vertical grid
    pad1[j]->Draw();             // Draw the Upper pad: pad1
    pad1[j]->cd();               // pad1 becomes the current pad
    pad1[j]->SetLogy();
    
    vFrame[j] = pad1[j]->DrawFrame(0, 0.3, 1.0, 700);
    histDraw(pad1[j],h1[j],h3[j],hs[j],vFrame[j]);
   
    h3[j]->Draw("ePsames");  //data 
    h2[j]->Draw("hist same");     //signal sample
    h9[j]->Draw("hist same");
    h6[j]->Draw("hist same");
    h7[j]->Draw("hist same");
    h8[j]->Draw("hist same");
    h10[j]->Draw("hist same");
    //h13[j]->Draw("sames");
    
    
    
    h[j]->Sumw2();
    h[j]->Add(h12[j]);
    h[j]->Add(h5[j]);
    h[j]->Add(h13[j]);
    h[j]->Add(h17[j]);
    h[j]->Add(h14[j]);
    h[j]->Add(h15[j]);
    h[j]->Add(h18[j]);
    h[j]->Add(h11[j]);
    h[j]->Add(h16[j]);
    h[j]->Add(h4[j]);
    h[j]->Add(h23[j]);
    h[j]->Add(h20[j]);
    h[j]->Add(h19[j]);
    h[j]->Add(h21[j]);
    h[j]->Add(h22[j]);
    sprintf(name,"h_err%i",j);
    h_err[j] = (TH1F*)h[j]->Clone(name);
    Int_t n=h_err[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++){
      if(h_err[j]->GetBinContent(i)!=0)
	h_err[j]->SetBinContent(i,1.0);
    }
    
    //h[j]->SetFillColorAlpha(40,0.35);
    //h[j]->SetMarkerStyle(20);
    h[j]->SetFillColor(12);
    //h[j]->SetLineWidth(20);
    h[j]->SetFillStyle(3018);
    h[j]->Draw("e2 same");
    
    vFrame[j]->GetYaxis()->SetTitle(ytitle);
    vFrame[j]->GetYaxis()->SetTitleSize(0.06);
    vFrame[j]->GetYaxis()->SetTitleOffset(0.8);
    vFrame[j]->GetYaxis()->SetLabelSize(0.05);
    vFrame[j]->GetXaxis()->SetLabelOffset(999);
    vFrame[j]->GetXaxis()->SetLabelSize(0);
    lg1->Draw();
    CMS_lumi( pad1[j], iPeriod, iPos );
    gPad->RedrawAxis();

    
    c[j]->cd();          // Go back to the main canvas before defining pad2

    sprintf(name,"pad2%i",j);
    pad2[j] = new TPad(name, name , 0, 0.0, 1.0, 0.3);
    pad2[j]->SetTopMargin(0.03);
    pad2[j]->SetBottomMargin(0.25);
    pad2[j]->SetGridx(); // vertical grid
    pad2[j]->SetGridy(); 
    pad2[j]->Draw();
    pad2[j]->cd();

   

    sprintf(name,"hd%i",j);
    hd[j] = (TH1F*)h3[j]->Clone(name); // For data
    
    hd[j]->SetLineColor(kBlack);
    hd[j]->Sumw2();
    hd[j]->SetStats(0);// No statistics on lower plot
    hd[j]->SetMinimum(0.5);  // Define Y ..
    hd[j]->SetMaximum(1.5); // .. range
    hd[j]->Divide(h[j]); //Data/MC
    hd[j]->SetMarkerStyle(21);
    hd[j]->Draw("ep");
    
    n=h[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      if(h[j]->GetBinContent(i)!=0){
	bin_err =h[j]->GetBinError(i)/h[j]->GetBinContent(i);
	h_err[j]->SetBinError(i,bin_err);
      }
    }
    h_err[j]->SetFillColor(40);
    //h[j]->SetLineWidth(20);
    h_err[j]->SetFillStyle(3001);
    h_err[j]->Draw("e2 same");
    hd[j]->Draw("ep same");

    
    hd[j]->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    
    hd[j]->GetYaxis()->SetTitle("Data/MC");
    hd[j]->GetYaxis()->SetNdivisions(310);
    hd[j]->GetYaxis()->SetTitleSize(30);
    hd[j]->GetYaxis()->SetTitleFont(43);
    hd[j]->GetYaxis()->SetTitleOffset(1.);
    hd[j]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[j]->GetYaxis()->SetLabelSize(20);

    // X axis ratio plot settings
    hd[j]->GetXaxis()->SetTitleSize(30);
    hd[j]->GetXaxis()->SetTitleFont(43);
    hd[j]->GetXaxis()->SetTitleOffset(2.8);
    hd[j]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hd[j]->GetXaxis()->SetLabelSize(25);
   
   sprintf(name,"l%i",j);
   
  }
  vFrame[0]->SetTitle("Dijet  P_{T} ");
  hd[0]->GetXaxis()->SetTitle("jj  P_{T} (GeV)");
  

  vFrame[1]->SetTitle("Dijet #eta ");
  hd[1]->GetXaxis()->SetTitle("jj #eta ");
  

  vFrame[2]->SetTitle("Dijet #phi ");
  hd[2]->GetXaxis()->SetTitle("jj #phi ");
  
  vFrame[3]->SetTitle("#Delta #eta between two leading leptons");
  hd[3]->GetXaxis()->SetTitle("#Delta #eta between #mu_{1} and #mu_{2}");
 
  vFrame[4]->SetTitle("#Delta R between two leading leptons ");
  hd[4]->GetXaxis()->SetTitle("#Delta R between #mu_{1} and #mu_{2} ");
  
  vFrame[5]->SetTitle("#Delta #phi between two leading leptons ");
  hd[5]->GetXaxis()->SetTitle("#Delta #phi between #mu_{1} and #mu_{2}");

  vFrame[6]->SetTitle("dilepton  P_{T} ");
  hd[6]->GetXaxis()->SetTitle("#mu#mu P_{T} (GeV) ");
  
  vFrame[7]->SetTitle("dilepton #eta ");
  hd[7]->GetXaxis()->SetTitle("#mu#mu #eta " );
  
  vFrame[8]->SetTitle("dilepton #phi ");
  hd[8]->GetXaxis()->SetTitle("#mu#mu #phi " );
  
  vFrame[9]->SetTitle("dilepton mass");
  hd[9]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");

  vFrame[10]->SetTitle("dilepton mass");
  hd[10]->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");

  vFrame[11]->SetTitle("Dijet mass");
  hd[11]->GetXaxis()->SetTitle("M(jj) (GeV)");

  vFrame[12]->SetTitle("#Delta #eta between j_{1} and j_{2}");
  hd[12]->GetXaxis()->SetTitle("#Delta #eta between j_{1} and j_{2}");

  vFrame[13]->SetTitle("#Delta #phi between j_{1} and j_{2}");
  hd[13]->GetXaxis()->SetTitle("#Delta #phi between j_{1} and j_{2}");

  vFrame[14]->SetTitle(" Leading Jet P_{T} ");
  hd[14]->GetXaxis()->SetTitle("j_{1} P_{T} (GeV) ");
 
  vFrame[15]->SetTitle("Leading Jet #eta ");
  hd[15]->GetXaxis()->SetTitle("j_{1}  #eta");
  
  vFrame[16]->SetTitle(" Sub-leading Jet P_{T} ");
  hd[16]->GetXaxis()->SetTitle("j_{2} P_{T} (GeV) ");
 
  vFrame[17]->SetTitle("Sub-leading Jet #eta ");
  hd[17]->GetXaxis()->SetTitle("j_{2}  #eta");

  vFrame[18]->SetTitle("cos(#theta^{*})");
  hd[18]->GetXaxis()->SetTitle("cos(#theta^{*})");

  vFrame[19]->SetTitle("Min #Delta R between muon and jet ");
  hd[19]->GetXaxis()->SetTitle("Min #Delta R (#mu j)");

  vFrame[20]->SetTitle("Max #Delta R between muon and jet ");
  hd[20]->GetXaxis()->SetTitle("Max #Delta R (#mu j)");

  vFrame[21]->SetTitle("Min #Delta R between dimuon and jet ");
  hd[21]->GetXaxis()->SetTitle("Min #Delta R (#mu#mu j) ");

  vFrame[22]->SetTitle("Max #Delta R between dimuon and jet ");
  hd[22]->GetXaxis()->SetTitle("Max #Delta R (#mu#mu j)");

  vFrame[23]->SetTitle("Zeppenfeld variable");
  hd[23]->GetXaxis()->SetTitle("Zeppenfeld variable");
  
  vFrame[24]->SetTitle("Dijet+dimuon mass");
  hd[24]->GetXaxis()->SetTitle("M(#mu#mu+jj) (GeV)");

  vFrame[25]->SetTitle(" Dimuon+Dijet P_{T} ");
  hd[25]->GetXaxis()->SetTitle("M(#mu#mu+jj) (GeV)");
  
  vFrame[26]->SetTitle("#Delta #eta of dimuon+dijet");
  hd[26]->GetXaxis()->SetTitle("#Delta #eta (#mu#mu+jj) ");

  vFrame[27]->SetTitle("#Delta #phi of dimuon+dijet");
  hd[27]->GetXaxis()->SetTitle("#Delta #phi (#mu#mu+jj)");

  vFrame[28]->SetTitle(" Leading Jet QGL ");
  hd[28]->GetXaxis()->SetTitle("jet 1 qgl");

  vFrame[29]->SetTitle(" sub-leading Jet QGL ");
  hd[29]->GetXaxis()->SetTitle("jet 2 qgl");
  vFrame[30]->SetTitle("# of EWK jets with P_{T} > 5 GeV ");
  hd[30]->GetXaxis()->SetTitle("# of soft EWK jets with P_{T} > 5 GeV");
  for(int j=0;j<31;j++){
    
    pad2[j]->cd();
    l[j]= new TLine(hd[j]->GetXaxis()->GetXmin(),1.0,hd[j]->GetXaxis()->GetXmax(),1.0);
    l[j]->SetLineColor(kBlack); l[j]->SetLineWidth(2); l[j]->SetLineStyle(1);
    l[j]->Draw("same");
   
   
    //c[j]->cd();
    //c[j]->Update();
    //c[j]->RedrawAxis();
    //c[j]->GetFrame()->Draw();
    /* name_unc=plotname[j]+".pdf";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".png";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".C";
    c[j]->SaveAs(name_unc);
    */

  }
  for(int j=0;j<31;j++){
    //c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".png";
    c[j]->SaveAs(name_unc);
    //name_unc=plotname[j]+".pdf";
    //c[j]->SaveAs(name_unc);
    //name_unc=plotname[j]+".C";
    //c[j]->SaveAs(name_unc);
  }
  

}


//Define the function decorate for histograms
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);
  
  h->SetMarkerSize(1.3);
  h->SetTitle(title);
  h->SetMinimum(0.00001);
}
//Define the function decorate for legends
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(100);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //g->SetHeader(legendheader);
}
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame )
{
  
  double maxCont = -1.0;
  maxCont=h1a->GetBinContent(h1a->GetMaximumBin());
  if(maxCont<h1b->GetBinContent(h1b->GetMaximumBin()))maxCont=h1b->GetBinContent(h1b->GetMaximumBin());
  maxCont*=100.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.01, maxRange, maxCont);
  hs->Draw("same hist");
}

// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
/*
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);

  return nevents;
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);

  return nevents;
}
*/
void h12ascii (TH1* h)
{
   Int_t n = h->GetNbinsX();
   
   for (Int_t i=1; i<=n; i++) {
      printf("%g %g\n",
             h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
             h->GetBinContent(i));
   }
}
