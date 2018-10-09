#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) ;
void histDraw(TPad *p,TH1F *h1a,TH1F *h1b, THStack *hs,TH1F *vFrame );

void DatavsMC()
{

  //First declare the file names
  TString f[77];
  f[0] = "DYJetsToLL.root"; //Fill in the appropriate names.
  f[1] = "VBFHToMuMu.root";
  f[2] = "Data2017.root";
  f[3] = "ttTosemileptonic.root";
  f[4] = "ttTo2l2v.root";
  f[5] = "ggH.root";
  f[6] = "WplusH.root";
  f[7] = "WminusH.root";
  f[8] = "ZH.root";
  f[9] = "ttH.root";
  f[10] = "WZTo1L1Nu2Q.root";
  f[11]="ZZ.root";
  f[12] = "WZTo3LNu.root";
  f[13]="WZTo2L2Q.root";
  f[14]="WWTo2L2Nu.root";
  f[15]="WWToLNuQQ.root";
  f[16] = "WWW_4F.root";
  f[17] = "WWZ_4F.root";
  /*f[18] = "ttbar_JECup_Ht500.root";
  f[19] = "TTPrime_800_JECup_Ht500_tZ.root";
  f[20] = "TTPrime_900_JECup_Ht500_tZ.root";
  f[21] = "TTPrime_1000_JECup_Ht500_tZ.root";
  f[22] = "TTPrime_1200_JECup_Ht500_tZ.root";
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
  plotname[0]= "mu1pt";
  plotname[1] = "mu1eta";
  plotname[2] = "mu1phi";
  plotname[3]= "mu2pt";
  plotname[4] = "mu2eta";
  plotname[5] = "mu2phi";
  plotname[6] ="mu1mu2dR";
  plotname[7] = "mu1mu2dPhi";
  plotname[8]= "diMuon_pt";
  plotname[9] = "diMuon_eta";
  plotname[10] = "diMuon_phi";
  plotname[11] = "diMuon_mass_SR";
  plotname[12] = "diMuon_mass_110To120";
  plotname[13]= "diMuon_mass_130To150";
  plotname[14] = "j1pt";
  plotname[15] = "j1phi";
  plotname[16] = "j1eta";
  plotname[17] = "j2pt";
  plotname[18] = "j2eta";
  plotname[19] = "j2phi";
  plotname[20] = "Njet";
  plotname[21] = "j1j2dR";
  plotname[22] = "j1j2dPhi";
  plotname[23] = "dijet_pt";
  plotname[24] = "dijet_eta";
  plotname[25] = "dijet_phi";
  plotname[26] = "diJet_mass";
  plotname[27] = "MET_pt";
  plotname[28] = "MET_phi";
  plotname[29] = "MET_sumEt";
  plotname[30]="diMuon_mass_110To150"; 
  plotname[31] = "Nbjet";
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
  for(int i=0;i<18;i++) file[i] = new TFile(f[i]);
  
  //change directory if required.
  //TDirectory *directory[77];
  //for(int i=0;i<77;i++) directory[i]= (TDirectory*) file[i]->Get("DatavsMC");


    /* Scaling can be done either by some outside measure (such as size of samples)
     or by normalizing the histograms such that they have the same integral.
     Conventionally, when comparing shapes, we shall normalize to 1
  */
    /*
  //==========for 30 fb-1======================//
  double scale_TT800=0.0657100711;
  double scale_TT900=0.029244122;//100% BR to tZtZ
  double scale_TT1000=0.0142894258;
  double scale_TT1100=0.0073985181;
  double scale_TT1200=0.0038505047;
  double scale_TT1300=0.0020999704;
  double scale_TT1400=0.001150546;
  double scale_ttbar=0.1370105511;
  double scale_DY_Ht100to200=0.6448872882;
  double scale_DY_Ht200to400=0.174180095;
  double scale_DY_Ht400to600=0.0251907386;
  //double scale_DY_Ht600toInf=0.0639509341;
  double scale_DY_Ht600to800=0.0061095439;
  double scale_DY_Ht800to1200=0.009408837;
  double scale_DY_Ht1200to2500=0.0069418046;
  double scale_DY_Ht2500toInf=0.0002541987;
    */
 //===========for 35.9 fb-1 data ===========//
 
 // double scale_ttbar=0.3866429988;
 



  double scale_DYJetsToLL= 5765.4*41.529*1000/(489144902631.884766+3241270753030.957031);
  double scale_VBFHToMuMu=0.000823*41.529*1000/4506449.599577;
  double scale_ttTosemileptonic=687.1*41.529*1000/10894741368.0;
  double scale_ttTo2l2v=85.656*41.529*1000/623402174.0;
  double scale_ggH=0.009605*1000*41.529/217554238.5;
  double scale_WplusH=0.000183*41.529*1000/259992.317749;
  double scale_WminusH=0.000116*41.529*1000/162196.811523;
  double scale_ZH=0.000192*41.529*1000/13507.2419434;
  double scale_ttH=0.000110*41.529*1000/155014.531738;
  double scale_WZTo1L1Nu2Q=11.61*41.529*1000/352741934.219;
  double scale_ZZ=16.523*41.529*1000/1949768.0;
  double scale_WZTo3LNu=4.42965*41.529*1000/93694769.25;
  double scale_WZTo2L2Q=5.595*41.529*1000/255256973.25;
  double scale_WWTo2L2Nu=12.46*41.529*1000/177178.179688;
  double scale_WWToLNuQQ=45.99*41.529*1000/405648754.016;
  double scale_WWW_4F=0.2086*41.529*1000/50039.244873;
  double scale_WWZ_4F=0.1651*41.529*1000/41205.3044434;

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
  TH1F *h[32], *h_err[32];
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
  
  for(int j=0;j<32;j++){
    h1[j] = (TH1F*)file[0]->Get(plotname[j]);
    h1[j]->Scale(scale_DYJetsToLL);
    /*
    name_unc=plotname[j]+"_Up";
    h1_Up[j] =(TH1F*)directory[0]->Get(name_unc);
    h1_Up[j]->Scale(scale_DY_Ht100to250);
    name_unc=plotname[j]+"_Down";
    h1_Down[j] =(TH1F*)directory[0]->Get(name_unc);
    h1_Down[j]->Scale(scale_DY_Ht100to250);
    */
    //===========================JEC Up
    /*h1_JECUp[j] = (TH1F*)directory[17]->Get(plotname[j]);
    h1_JECUp[j]->Scale(scale_DY_Ht100to250);

    h1_JECDown[j] = (TH1F*)directory[32]->Get(plotname[j]);
    h1_JECDown[j]->Scale(scale_DY_Ht100to250);


    h1_JERUp[j] = (TH1F*)directory[47]->Get(plotname[j]);
    h1_JERUp[j]->Scale(scale_DY_Ht100to250);

    h1_JERDown[j] = (TH1F*)directory[62]->Get(plotname[j]);
    h1_JERDown[j]->Scale(scale_DY_Ht100to250);
    
    Int_t n=h1_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      //cout<<j<<endl;
      bin_err =(h1_Up[j]->GetBinContent(i)-h1_Down[j]->GetBinContent(i));
      jec_err =(h1_JECUp[j]->GetBinContent(i)-h1_JECDown[j]->GetBinContent(i));
      jer_err =(h1_JERUp[j]->GetBinContent(i)-h1_JERDown[j]->GetBinContent(i));*/
    //stat_err = h1[j]->GetBinError(i);
      //cout<<bin_err<<" "<<stat_err<<endl;
      //bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      //cout<<plotname[j]<<" "<<i<<" "<<bin_err<<endl;
      //h1[j]->SetBinError(i,bin_err);
    //}

    sprintf(name,"h%i",j);
    h[j] = (TH1F*)h1[j]->Clone(name); //For MC // doing it here so that I can adjust the color of h[j]
    decorate(h1[j],"",ytitle,"",kMagenta-9,2,kMagenta-9,20,1);
   if(j==11)cout<<h1[j]->Integral()<<endl;
  }

  for(int j=0;j<32;j++){
    h2[j] = (TH1F*)file[1]->Get(plotname[j]);
    /*
    name_unc=plotname[j]+"_Up";
    h2_Up[j] =(TH1F*)directory[1]->Get(name_unc);
    h2_Up[j]->Scale(scale_ttbar);
    name_unc=plotname[j]+"_Down";
    h2_Down[j] =(TH1F*)directory[1]->Get(name_unc);
    h2_Down[j]->Scale(scale_ttbar);

    h2_JECUp[j] = (TH1F*)directory[18]->Get(plotname[j]);
    h2_JECUp[j]->Scale(scale_ttbar);

    h2_JECDown[j] = (TH1F*)directory[33]->Get(plotname[j]);
    h2_JECDown[j]->Scale(scale_ttbar);

    h2_JERUp[j] = (TH1F*)directory[48]->Get(plotname[j]);
    h2_JERUp[j]->Scale(scale_ttbar);

    h2_JERDown[j] = (TH1F*)directory[63]->Get(plotname[j]);
    h2_JERDown[j]->Scale(scale_ttbar);
    Int_t n=h2_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h2_Up[j]->GetBinContent(i)-h2_Down[j]->GetBinContent(i));
      jec_err =(h2_JECUp[j]->GetBinContent(i)-h2_JECDown[j]->GetBinContent(i));
      jer_err =(h2_JERUp[j]->GetBinContent(i)-h2_JERDown[j]->GetBinContent(i));
      stat_err = h2[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h2[j]->SetBinError(i,bin_err);
      }*/
     decorate(h2[j],"",ytitle,"",kGreen-9,2,kGreen-9,21,0);
     h2[j]->Scale(scale_VBFHToMuMu);
    if(j==11)cout<<h2[j]->Integral()<<endl;
  }

  for(int j=0;j<32;j++){
    h3[j] = (TH1F*)file[2]->Get(plotname[j]);
    decorate(h3[j],"",ytitle,"",kBlack,2,kBlack,8,0);
  }
  
  for(int j=0;j<32;j++){
    h4[j] = (TH1F*)file[3]->Get(plotname[j]);
    decorate(h4[j],"",ytitle,"",kGreen+2,2,kGreen+2,21,1);

    h4[j]->Scale(scale_ttTosemileptonic);
    if(j==11)cout<<h4[j]->Integral()<<endl;
  }
 
  for(int j=0;j<32;j++){
    h5[j] = (TH1F*)file[4]->Get(plotname[j]);
    //cout<<"tt 2l2v\n";
    decorate(h5[j],"",ytitle,"",kYellow+2,2,kYellow+2,21,1);
    h5[j]->Scale(scale_ttTo2l2v);
    if(j==11)cout<<h5[j]->Integral()<<endl;
  }
  cout<<"scale ggh:"<<scale_ggH<<endl;
  for(int j=0;j<32;j++){
    h6[j] = (TH1F*)file[5]->Get(plotname[j]);
    //cout<<"ggH\n";
    decorate(h6[j],"",ytitle,"",kCyan+4,2,kCyan+4,21,0);
    h6[j]->Scale(scale_ggH);
    
    if(j==11){
       cout<<h6[j]->Integral()<<endl;
    }
    
  }

  for(int j=0;j<32;j++){
    h7[j] = (TH1F*)file[6]->Get(plotname[j]);
    decorate(h7[j],"",ytitle,"",kBlue,2,kBlue,21,0);
    h7[j]->Scale(scale_WplusH);
    if(j==11)cout<<h7[j]->Integral()<<endl;
  }

  for(int j=0;j<32;j++){
    h8[j] = (TH1F*)file[7]->Get(plotname[j]);
    //cout<<h8[j]->Integral();
    decorate(h8[j],"",ytitle,"",kPink-7,2,kPink-7,21,0);
    h8[j]->Scale(scale_WminusH);
    if(j==11)cout<<h8[j]->Integral()<<endl;
    /*
    name_unc=plotname[j]+"_Up";
    h8_Up[j] =(TH1F*)directory[7]->Get(name_unc);
    h8_Up[j]->Scale(scale_DY_Ht250to400);
    name_unc=plotname[j]+"_Down";
    h8_Down[j] =(TH1F*)directory[7]->Get(name_unc);
    h8_Down[j]->Scale(scale_DY_Ht250to400);

    h8_JECUp[j] = (TH1F*)directory[23]->Get(plotname[j]);
    h8_JECUp[j]->Scale(scale_DY_Ht250to400);
    
    h8_JECDown[j] = (TH1F*)directory[38]->Get(plotname[j]);
    h8_JECDown[j]->Scale(scale_DY_Ht250to400);

    h8_JERUp[j] = (TH1F*)directory[53]->Get(plotname[j]);
    h8_JERUp[j]->Scale(scale_DY_Ht250to400);

    h8_JERDown[j] = (TH1F*)directory[68]->Get(plotname[j]);
    h8_JERDown[j]->Scale(scale_DY_Ht250to400);
    Int_t n=h1_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h8_Up[j]->GetBinContent(i)-h8_Down[j]->GetBinContent(i));
      stat_err = h8[j]->GetBinError(i);
      jec_err =(h8_JECUp[j]->GetBinContent(i)-h8_JECDown[j]->GetBinContent(i));
      jer_err =(h8_JERUp[j]->GetBinContent(i)-h8_JERDown[j]->GetBinContent(i));
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h8[j]->SetBinError(i,bin_err);
      }*/
  }
    
  
  for(int j=0;j<32;j++){
    h9[j] = (TH1F*)file[8]->Get(plotname[j]);
    decorate(h9[j],"",ytitle,"",kRed-2,2,kRed-2,21,0);
    //cout<<"Zh\n";
    h9[j]->Scale(scale_ZH);
    if(j==11)cout<<h9[j]->Integral()<<endl;
    /*  
    name_unc=plotname[j]+"_Up";
    h9_Up[j] =(TH1F*)directory[8]->Get(name_unc);
    h9_Up[j]->Scale(scale_DY_Ht400to650);
    name_unc=plotname[j]+"_Down";
    h9_Down[j] =(TH1F*)directory[8]->Get(name_unc);
    h9_Down[j]->Scale(scale_DY_Ht400to650);

    h9_JECUp[j] = (TH1F*)directory[24]->Get(plotname[j]);
    h9_JECUp[j]->Scale(scale_DY_Ht400to650);
    
    h9_JECDown[j] = (TH1F*)directory[39]->Get(plotname[j]);
    h9_JECDown[j]->Scale(scale_DY_Ht400to650);
    
    h9_JERUp[j] = (TH1F*)directory[54]->Get(plotname[j]);
    h9_JERUp[j]->Scale(scale_DY_Ht400to650);

    h9_JERDown[j] = (TH1F*)directory[69]->Get(plotname[j]);
    h9_JERDown[j]->Scale(scale_DY_Ht400to650);
    Int_t n=h9_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h9_Up[j]->GetBinContent(i)-h9_Down[j]->GetBinContent(i));
      jec_err =(h9_JECUp[j]->GetBinContent(i)-h9_JECDown[j]->GetBinContent(i));
      jer_err =(h9_JERUp[j]->GetBinContent(i)-h9_JERDown[j]->GetBinContent(i));
      stat_err = h9[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h9[j]->SetBinError(i,bin_err);
      }*/
      }
  
  for(int j=0;j<32;j++){
    h10[j] = (TH1F*)file[9]->Get(plotname[j]);
    decorate(h10[j],"",ytitle,"",kRed,2,kRed,21,0);
    h10[j]->Scale(scale_ttH);
    if(j==11)cout<<h10[j]->Integral()<<endl;
    /*
    name_unc=plotname[j]+"_Up";
    h10_Up[j] =(TH1F*)directory[9]->Get(name_unc);
    h10_Up[j]->Scale(scale_DY_Ht650toInf);
    name_unc=plotname[j]+"_Down";
    h10_Down[j] =(TH1F*)directory[9]->Get(name_unc);
    h10_Down[j]->Scale(scale_DY_Ht650toInf);

    h10_JECUp[j] = (TH1F*)directory[25]->Get(plotname[j]);
    h10_JECUp[j]->Scale(scale_DY_Ht650toInf);
    
    h10_JECDown[j] = (TH1F*)directory[40]->Get(plotname[j]);
    h10_JECDown[j]->Scale(scale_DY_Ht650toInf);
    
    h10_JERUp[j] = (TH1F*)directory[55]->Get(plotname[j]);
    h10_JERUp[j]->Scale(scale_DY_Ht650toInf);

    h10_JERDown[j] = (TH1F*)directory[70]->Get(plotname[j]);
    h10_JERDown[j]->Scale(scale_DY_Ht650toInf);
    Int_t n=h10_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h10_Up[j]->GetBinContent(i)-h10_Down[j]->GetBinContent(i));
      jec_err =(h10_JECUp[j]->GetBinContent(i)-h10_JECDown[j]->GetBinContent(i));
      jer_err =(h10_JERUp[j]->GetBinContent(i)-h10_JERDown[j]->GetBinContent(i));
      stat_err = h10[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h10[j]->SetBinError(i,bin_err);
      }*/
  }

  
  for(int j=0;j<32;j++){
    h11[j] = (TH1F*)file[10]->Get(plotname[j]);
    decorate(h11[j],"",ytitle,"",kOrange-8,2,kOrange-8,21,1);
    h11[j]->Scale(scale_WZTo1L1Nu2Q);
    if(j==11)cout<<h11[j]->Integral()<<endl;
  }
    for(int j=0;j<32;j++){
    h12[j] = (TH1F*)file[11]->Get(plotname[j]);
    decorate(h12[j],"",ytitle,"",kYellow-9,2,kYellow-9,21,1);
    h12[j]->Scale(scale_ZZ);
    if(j==11)cout<<h12[j]->Integral()<<endl;
  }
  
  for(int j=0;j<32;j++){
    h13[j] = (TH1F*)file[12]->Get(plotname[j]);
    decorate(h13[j],"",ytitle,"",kBlue-10,2,kBlue-10,21,1);
    h13[j]->Scale(scale_WZTo3LNu);
    if(j==11)cout<<h13[j]->Integral()<<endl;
  }

for(int j=0;j<32;j++){
    h14[j] = (TH1F*)file[13]->Get(plotname[j]);
    decorate(h14[j],"",ytitle,"",kRed-7,2,kRed-7,21,1);
    h14[j]->Scale(scale_WZTo2L2Q);
    if(j==11)cout<<h14[j]->Integral()<<endl;
    /*
    name_unc=plotname[j]+"_Up";
    h14_Up[j] =(TH1F*)directory[13]->Get(name_unc);
    h14_Up[j]->Scale(scale_WW);
    name_unc=plotname[j]+"_Down";
    h14_Down[j] =(TH1F*)directory[13]->Get(name_unc);
    h14_Down[j]->Scale(scale_WW);

    h14_JECUp[j] = (TH1F*)directory[29]->Get(plotname[j]);
    h14_JECUp[j]->Scale(scale_WW);

    h14_JECDown[j] = (TH1F*)directory[44]->Get(plotname[j]);
    h14_JECDown[j]->Scale(scale_WW);

    h14_JERUp[j] = (TH1F*)directory[59]->Get(plotname[j]);
    h14_JERUp[j]->Scale(scale_WW);

    h14_JERDown[j] = (TH1F*)directory[74]->Get(plotname[j]);
    h14_JERDown[j]->Scale(scale_WW);
    
    Int_t n=h14_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h14_Up[j]->GetBinContent(i)-h14_Down[j]->GetBinContent(i));
      jec_err =(h14_JECUp[j]->GetBinContent(i)-h14_JECDown[j]->GetBinContent(i));
      jer_err =(h14_JERUp[j]->GetBinContent(i)-h14_JERDown[j]->GetBinContent(i));
      stat_err = h14[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h14[j]->SetBinError(i,bin_err);
      }*/
  }
  for(int j=0;j<32;j++){
    h15[j] = (TH1F*)file[14]->Get(plotname[j]);
    decorate(h15[j],"",ytitle,"",kBlue-6,2,kBlue-6,21,1);
    h15[j]->Scale(scale_WWTo2L2Nu);
    if(j==11)cout<<h15[j]->Integral()<<endl;
    /*
    name_unc=plotname[j]+"_Up";
    h15_Up[j] =(TH1F*)directory[14]->Get(name_unc);
    h15_Up[j]->Scale(scale_WZ);
    name_unc=plotname[j]+"_Down";
    h15_Down[j] =(TH1F*)directory[14]->Get(name_unc);
    h15_Down[j]->Scale(scale_WZ);

    h15_JECUp[j] = (TH1F*)directory[30]->Get(plotname[j]);
    h15_JECUp[j]->Scale(scale_WZ);

    h15_JECDown[j] = (TH1F*)directory[45]->Get(plotname[j]);
    h15_JECDown[j]->Scale(scale_WZ);

    h15_JERUp[j] = (TH1F*)directory[60]->Get(plotname[j]);
    h15_JERUp[j]->Scale(scale_WZ);

    h15_JERDown[j] = (TH1F*)directory[75]->Get(plotname[j]);
    h15_JERDown[j]->Scale(scale_WZ);
    Int_t n=h15_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h15_Up[j]->GetBinContent(i)-h15_Down[j]->GetBinContent(i));
      jec_err =(h15_JECUp[j]->GetBinContent(i)-h15_JECDown[j]->GetBinContent(i));
      jer_err =(h15_JERUp[j]->GetBinContent(i)-h15_JERDown[j]->GetBinContent(i));
      stat_err = h15[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h15[j]->SetBinError(i,bin_err);
      }*/
  }
  
  for(int j=0;j<32;j++){
    h16[j] = (TH1F*)file[15]->Get(plotname[j]);
    decorate(h16[j],"",ytitle,"",kViolet+1,2,kViolet+1,21,1);
    h16[j]->Scale(scale_WWToLNuQQ);
    if(j==11)cout<<h16[j]->Integral()<<endl;
    /*  
    name_unc=plotname[j]+"_Up";
    h16_Up[j] =(TH1F*)directory[15]->Get(name_unc);
    h16_Up[j]->Scale(scale_WWW_4F);
    name_unc=plotname[j]+"_Down";
    h16_Down[j] =(TH1F*)directory[15]->Get(name_unc);
    h16_Down[j]->Scale(scale_ZZ);

    h16_JECUp[j] = (TH1F*)directory[31]->Get(plotname[j]);
    h16_JECUp[j]->Scale(scale_ZZ);

    h16_JECDown[j] = (TH1F*)directory[46]->Get(plotname[j]);
    h16_JECDown[j]->Scale(scale_ZZ);

    h16_JERUp[j] = (TH1F*)directory[61]->Get(plotname[j]);
    h16_JERUp[j]->Scale(scale_ZZ);

    h16_JERDown[j] = (TH1F*)directory[76]->Get(plotname[j]);
    h16_JERDown[j]->Scale(scale_ZZ);
    Int_t n=h16_Up[j]->GetNbinsX();
    for (Int_t i=1; i<=n; i++) {
      bin_err =(h16_Up[j]->GetBinContent(i)-h16_Down[j]->GetBinContent(i));
      jec_err =(h16_JECUp[j]->GetBinContent(i)-h16_JECDown[j]->GetBinContent(i));
      jer_err =(h16_JERUp[j]->GetBinContent(i)-h16_JERDown[j]->GetBinContent(i));
      stat_err = h16[j]->GetBinError(i);
      bin_err = sqrt(bin_err*bin_err+stat_err*stat_err+jec_err*jec_err+jer_err*jer_err+0.026*0.026+0.03*0.03+0.15*0.15+0.01*0.01);
      h16[j]->SetBinError(i,bin_err);
      }*/
  }
  for(int j=0;j<32;j++){
    h17[j] = (TH1F*)file[16]->Get(plotname[j]);
    decorate(h17[j],"",ytitle,"",kRed-10,2,kRed-10,21,1);
    h17[j]->Scale(scale_WWW_4F);
    if(j==11)cout<<h17[j]->Integral()<<endl;
  }

  for(int j=0;j<32;j++){
    h18[j] = (TH1F*)file[17]->Get(plotname[j]);
    decorate(h18[j],"",ytitle,"",kGray+2,2,kGray+2,21,1);
    h18[j]->Scale(scale_WWZ_4F);
    if(j==11)cout<<h18[j]->Integral()<<endl;
  }

  
  TLegend *lg1 = new TLegend(0.905,0.25,0.95,0.85);
  lg1->SetTextFont(132);
  lg1->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.02);
  lg1->AddEntry((TObject*)0,"Inclusive","");
  lg1->AddEntry(h3[3],"Data:","p");
  lg1->AddEntry((TObject*)0,"41.529 fb^{-1}",""); 
  
  lg1->AddEntry(h1[3],"DY+jets","f");
  lg1->AddEntry(h4[3],"t#bar{t} semi-leptonic","f");
  lg1->AddEntry(h5[3],"t#bar{t} #rightarrow 2l 2#nu","f");
  lg1->AddEntry(h11[3],"WZ #rightarrow 1l 1#nu 2Q","f");
  lg1->AddEntry(h13[3],"WZ #rightarrow 3l 1#nu","f");
  lg1->AddEntry(h14[3],"WZ #rightarrow 2l 2Q","f");
  lg1->AddEntry(h15[3],"WW #rightarrow 2l 2#nu","f");
  lg1->AddEntry(h16[3],"WW #rightarrow 1l 1#nu 2Q","f");
  lg1->AddEntry(h12[3],"ZZ","f");
  lg1->AddEntry(h17[3],"WWW #rightarrow 4F","f");
  lg1->AddEntry(h18[3],"WWZ #rightarrow 4F","f");
  
  lg1->AddEntry(h2[3],"VBF","l");
  lg1->AddEntry(h6[3],"ggH","l");
  lg1->AddEntry(h7[3],"W^{+}H","l");
  lg1->AddEntry(h8[3],"W^{-}H","l");
  lg1->AddEntry(h10[3],"ttH","l");
  lg1->AddEntry(h9[3],"ZH","l");
  for(int j=0;j<32;j++){
   
    sprintf(name,"hs%i",j);
    hs[j] = new THStack(name,plotname[j]); // Stack of all backgrounds
    hs[j]->Add(h18[j]);
    hs[j]->Add(h17[j]);
    hs[j]->Add(h13[j]);
    hs[j]->Add(h14[j]);
    hs[j]->Add(h11[j]);
    hs[j]->Add(h15[j]);
    hs[j]->Add(h12[j]);
    hs[j]->Add(h16[j]);

    hs[j]->Add(h5[j]);
    hs[j]->Add(h4[j]); //ttsl
   
    hs[j]->Add(h1[j]);
    

    sprintf(name,"c%i",j);
    c[j] = new TCanvas(name,name,800,800);
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
  vFrame[0]->SetTitle("Leading Lepton  P_{T} ");
  hd[0]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[1]->SetTitle("Leading Lepton #eta ");
  hd[1]->GetXaxis()->SetTitle(xtitle1);
  

  vFrame[2]->SetTitle("Leading Lepton #phi ");
  hd[2]->GetXaxis()->SetTitle(xtitle2);
  

  vFrame[3]->SetTitle("Sub-Leading Lepton  P_{T} ");
  hd[3]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[4]->SetTitle("Sub-Leading Lepton #eta ");
  hd[4]->GetXaxis()->SetTitle(xtitle1);
  

  vFrame[5]->SetTitle("Sub-Leading Lepton #phi ");
  hd[5]->GetXaxis()->SetTitle(xtitle2);
 
  vFrame[6]->SetTitle("#Delta R between two leading leptons ");
  hd[6]->GetXaxis()->SetTitle("#Delta R");
  

  vFrame[7]->SetTitle("#Delta #phi between two leading leptons ");
  hd[7]->GetXaxis()->SetTitle("#Delta #phi");

  vFrame[8]->SetTitle("dilepton  P_{T} ");
  hd[8]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[9]->SetTitle("dilepton #eta ");
  hd[9]->GetXaxis()->SetTitle(xtitle1);
  
  vFrame[10]->SetTitle("dilepton #phi ");
  hd[10]->GetXaxis()->SetTitle(xtitle2);
  
  vFrame[11]->SetTitle("dilepton mass");
  hd[11]->GetXaxis()->SetTitle("Mass (GeV)");

  vFrame[12]->SetTitle("dilepton mass");
  hd[12]->GetXaxis()->SetTitle("Mass (GeV)");

  vFrame[13]->SetTitle("dilepton mass");
  hd[13]->GetXaxis()->SetTitle("Mass (GeV)");


  
  vFrame[14]->SetTitle("Leading Jet  P_{T} ");
  hd[14]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[15]->SetTitle("Leading Jet #eta ");
  hd[15]->GetXaxis()->SetTitle(xtitle1);
  

  vFrame[16]->SetTitle("Leading Jet #phi ");
  hd[16]->GetXaxis()->SetTitle(xtitle2);
  

  vFrame[17]->SetTitle("Sub-Leading Jet  P_{T} ");
  hd[17]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[18]->SetTitle("Sub-Leading Jet #eta ");
  hd[18]->GetXaxis()->SetTitle(xtitle1);
  

  vFrame[19]->SetTitle("Sub-Leading Jet #phi ");
  hd[19]->GetXaxis()->SetTitle(xtitle2);

  vFrame[20]->SetTitle("Number of Jets");
  hd[20]->GetXaxis()->SetTitle("Njets");

  vFrame[21]->SetTitle("#Delta R between two leading jets");
  hd[21]->GetXaxis()->SetTitle("#Delta R");
  
  vFrame[22]->SetTitle("#Delta #phi between two leading jets ");
  hd[22]->GetXaxis()->SetTitle("#Delta #phi");

 
  vFrame[23]->SetTitle("Dijet P_{T} ");
  hd[23]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[24]->SetTitle("Dijet #eta ");
  hd[24]->GetXaxis()->SetTitle(xtitle1);
  
  vFrame[25]->SetTitle("Dijet #phi ");
  hd[25]->GetXaxis()->SetTitle(xtitle2);

  vFrame[26]->SetTitle("Dijet mass ");
  hd[26]->GetXaxis()->SetTitle("Mass (GeV)");

  vFrame[27]->SetTitle("MET P_{T} ");
  hd[27]->GetXaxis()->SetTitle(xtitle);
  

  vFrame[28]->SetTitle("MET #phi ");
  hd[28]->GetXaxis()->SetTitle(xtitle2);
  
  vFrame[29]->SetTitle("MET Sum Et ");
  hd[29]->GetXaxis()->SetTitle("Sum Et (GeV)");
  

  vFrame[30]->SetTitle("dilepton mass");
  hd[30]->GetXaxis()->SetTitle("Mass (GeV)");
  
  vFrame[31]->SetTitle("Number of b jets");
  hd[31]->GetXaxis()->SetTitle("Nbjet");
  
  for(int j=0;j<32;j++){
    c[j]->cd();
    c[j]->Update();
    pad2[j]->cd();
    l[j]= new TLine(hd[j]->GetXaxis()->GetXmin(),1.0,hd[j]->GetXaxis()->GetXmax(),1.0);
    l[j]->SetLineColor(kBlack); l[j]->SetLineWidth(2); l[j]->SetLineStyle(1);
    l[j]->Draw("same");

    name_unc=plotname[j]+".pdf";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".png";
    c[j]->SaveAs(name_unc);
    name_unc=plotname[j]+".C";
    c[j]->SaveAs(name_unc);


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
  g->SetFillColor(10);
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
  maxCont*=5.0;
  double maxRange = h1a->GetXaxis()->GetXmax();
  double minRange=h1a->GetXaxis()->GetXmin();
  //if(h1a->GetBinCenter(h1a->GetNbinsX()) > maxRange) maxRange =  h1a->GetBinCenter(h1a->GetNbinsX());
  vFrame = p->DrawFrame(minRange, 0.095, maxRange, maxCont);
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
