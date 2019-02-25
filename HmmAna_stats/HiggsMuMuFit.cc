#define HiggsMuMuFit_cxx
#include "HiggsMuMuFit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

using namespace RooFit ;

int main(int argc, char* argv[])
{

  if(argc < 3) {
    cerr << "Please give 4 arguments: input " << " " << " outputFileName" << " " << "is bkg or not " <<endl;
    return -1;
  }
  TString inputFileList = argv[1];
  TString outFileName = argv[2];
  TString isbkgstr = argv[3];

  TString mass_cut = "Higgs_mass > 110 && Higgs_mass< 150 &&";
  //e.g. "(cat_index==5)&&disc>0.97&&disc_ttbar>0.8&&(Higgs_mass>110&&Higgs_mass<150)";
  
  TString selection = mass_cut + " cat_index<4"; //ttH
  int category = 1;
  TString proc = isbkgstr;
  HiggsMuMuFit hmm(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf(); 

  /* 
  selection = "Higgs_mass > 110 && Higgs_mass< 150 && cat_index==2"; 
  category = 2;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

  selection = "Higgs_mass > 110 && Higgs_mass< 150 && cat_index==3"; 
  category = 3;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();
  */

  
  selection = mass_cut + " cat_index>3 && cat_index<7"; //VH-lep
  category = 2;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

  selection = mass_cut + " cat_index==7"; //VH-had? VBF?
  category = 3;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

  selection = mass_cut + " cat_index==8"; //VBF? VH-had??
  category = 4;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();
  
  selection = mass_cut + " cat_index>8 && fabs(reco_mu1_eta)<0.9 && disc_nw2016 >0.0";
  category = 5;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();
  
  selection = mass_cut + " cat_index>8 && fabs(reco_mu1_eta)<0.9 && disc_nw2016 <0.0";
  category = 6;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

  selection = mass_cut + " cat_index>8 && fabs(reco_mu1_eta)>0.9 && fabs(reco_mu1_eta)<1.9 && disc_nw2016>0.0";
  category = 7;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();


  selection = mass_cut + " cat_index>8 && fabs(reco_mu1_eta)>0.9 && fabs(reco_mu1_eta)<1.9 && disc_nw2016 <0.0";
  category = 8;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

 
  selection = mass_cut + " cat_index>8 && fabs(reco_mu1_eta)>1.9";
  category = 9;
  proc = isbkgstr;
  hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();
   

  return 0;
}

void HiggsMuMuFit::sigfit()
{ 
  cout <<"sigfit"<<endl;


  RooRealVar* catin= new RooRealVar("cat_index","cat_index",1,0,100) ;
  RooRealVar* reco_mu1_eta = new RooRealVar("reco_mu1_eta","reco_mu1_eta",1,-100,100) ;
  RooRealVar* disc_nw2016 = new RooRealVar("disc_nw2016","disc_nw2016",0,-1,1) ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",110,110,150) ;
  RooRealVar* cut_based_ct = new RooRealVar("cat_index","cat_index",1,0,10) ;
  RooRealVar* disc = new RooRealVar("disc","disc",0,0,1) ;
  RooRealVar* disc_ttbar = new RooRealVar("disc_ttbar","disc_ttbar",0,0,1) ;

  RooRealVar* mean1 = new RooRealVar(Form(procname+"cat%d_mean1",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex),125, 120, 130) ;
  RooRealVar* mean2 = new RooRealVar(Form(procname+"cat%d_mean2",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex),125, 120, 130) ;
  RooRealVar* mean3 = new RooRealVar(Form(procname+"cat%d_mean3",catindex),Form("mean3 of gaussians for "+procname+"cat%d", catindex),125, 120, 130) ;

  RooRealVar* sigma1 = new RooRealVar(Form(procname+"cat%d_sigma1",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex),1, 0.0, 100) ;
  RooRealVar* sigma2 = new RooRealVar(Form(procname+"cat%d_sigma2",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex),1, 0.0, 100) ;
  RooRealVar* sigma3 = new RooRealVar(Form(procname+"cat%d_sigma3",catindex),Form("sigma3 of gaussians for "+procname+"cat%d", catindex),1, 0.0, 100) ;

  RooGaussian* sig1 = new RooGaussian(Form(procname+"cat%d_sig1",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
  RooGaussian* sig2 = new RooGaussian(Form(procname+"cat%d_sig2",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;
  RooGaussian* sig3 = new RooGaussian(Form(procname+"cat%d_sig3",catindex),Form("Signal component 3 for "+procname+"cat%d", catindex),*mdimu, *mean3, *sigma3) ;

  RooRealVar* sig1frac = new RooRealVar(Form(procname+"cat%d_sig1frac",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),0.8,0.,1.) ;
  RooRealVar* sig2frac = new RooRealVar(Form(procname+"cat%d_sig2frac",catindex),Form("signal fraction of component 2 for "+procname+"cat%d", catindex),0.1,0.,1.) ;
  RooAddPdf* dimupdf = new RooAddPdf(Form(procname+"cat%d_sig_pdf",catindex),Form(procname+"cat%d_sig_pdf",catindex),RooArgList(*sig1,*sig2,*sig3),RooArgList(*sig1frac,*sig2frac)) ;

  RooRealVar* nsig = new RooRealVar(Form(procname+"cat%d_nsig",catindex),Form("number of signal events in signalRange "+procname+"cat%d", catindex),500,0,1000000) ;
  RooExtendPdf* sigpdf = new RooExtendPdf(Form(procname+"cat%d_sig_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*dimupdf,*nsig,"sigRange") ;

  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-100000,100000) ;

  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  obsAndWeight.add(*catin);
  obsAndWeight.add(*disc_nw2016);
  obsAndWeight.add(*reco_mu1_eta);
  obsAndWeight.add(*evWeight);

  cout <<"cut: "<<cut<<endl;

  string dataset_ch_name(Form("Sig"+procname+"_cat%d",catindex));
  RooCmdArg cut_arg = cut==""?RooCmdArg::none():RooFit::Cut(cut);
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(),
                RooArgSet(obsAndWeight), RooFit::WeightVar(*evWeight),
                RooFit::Import(*data_chain),
                cut_arg
                );

  RooFitResult* r_full = sigpdf->fitTo(*data,Save()) ;

  RooFitResult* r = sigpdf->fitTo(*data,RooFit::Strategy(1), RooFit::Extended(kTRUE), RooFit::Save(kTRUE));

  cout << "result of fit on all data " << endl ;
  //r_full->Print() ;
  cout << "result of fit in in signal region (note increased error on signal fraction)" << endl ;
  //r->Print() ;


  sig1frac->setConstant(true);
  sig2frac->setConstant(true);
  sigma1->setConstant(true);
  sigma2->setConstant(true);
  sigma3->setConstant(true);
  mean1->setConstant(true);
  mean2->setConstant(true);
  mean3->setConstant(true);
  nsig->setConstant(true);

  RooBinning tbins(110,150);
  tbins.addUniform(40,110,150) ;
  RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
  data->plotOn(dtframe,Binning(tbins));
  sigpdf->plotOn(dtframe);
  sigpdf->plotOn(dtframe, Components(*sig1), LineColor(kRed), LineStyle(kDashed));
  sigpdf->plotOn(dtframe, Components(*sig2), LineColor(kGreen), LineStyle(kDashed));
  sigpdf->plotOn(dtframe, Components(*sig3), LineColor(kOrange), LineStyle(kDashed));

  sigpdf->paramOn(dtframe,Layout(0.43, 0.93, 0.75), Format("NEU",AutoPrecision(0))) ;
  //dimupdf->paramOn(dtframe,Layout(0.43, 0.93, 0.65), Format("NEU",AutoPrecision(0))) ;
 
  data->statOn(dtframe,Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(0))) ;
  //Form(procname+"cat%d_sig_pdf",catindex)//
  /*
  dtframe->getAttText(Form(procname+"cat%d_sig_pdf",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_sig_pdf",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_sig_pdf",catindex))->SetTextSize(0.032);
  */
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextSize(0.032);
  dtframe->getAttText()->SetTextSize(0.032);

  TCanvas* c1 = new TCanvas();
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m.png",catindex)));

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;
  w->import(*sigpdf) ;
  w->import(*data) ;
  w->Print() ;
  w->writeToFile(TString::Format(outfile+"/pdfs/Hmm.input"+procname+"sig_13TeV.cat%d_m.root",catindex));

  delete data_chain;
  delete data;
  delete sigpdf; 
}   

void HiggsMuMuFit::bkgfit()
{
  cout <<"bkgfit"<<endl;    

  RooRealVar* catin= new RooRealVar("cat_index","cat_index",1,0,100) ;
  RooRealVar* reco_mu1_eta = new RooRealVar("reco_mu1_eta","reco_mu1_eta",1,-100,100) ;
  RooRealVar* disc_nw2016 = new RooRealVar("disc_nw2016","disc_nw2016",0,-1,1) ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",110,110,150) ;
  RooRealVar* cut_based_ct = new RooRealVar("cat_index","cat_index",1,0,10) ;
  RooRealVar* disc = new RooRealVar("disc","disc",0,0,1) ;
  RooRealVar* disc_ttbar = new RooRealVar("disc_ttbar","disc_ttbar",0,0,1) ;

  RooRealVar* mdimuslope0 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_slope0_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_slope0_cat%d",catindex),1,-50,1000) ;
  RooRealVar* mdimuslope1 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_slope1_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_slope1_cat%d",catindex),1,-50,1000) ;
  RooRealVar* mdimuslope2 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_slope2_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_slope2_cat%d",catindex),1,-50,1000) ;
  RooAbsPdf* mdimupdf = new RooExponential(TString::Format(procname+"cat%d_bkg_pdf",catindex),TString::Format(procname+"cat%d_bkg_pdf",catindex), *mdimu, *mdimuslope0);

  RooRealVar* nbkg = new RooRealVar(Form(procname+"cat%d_nbkg",catindex),Form("number of background events in signalRange "+procname+"cat%d", catindex),500,0,10000000) ;
  RooExtendPdf* bkgpdf = new RooExtendPdf(Form(procname+"cat%d_bkg_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*mdimupdf,*nbkg,"bkgRange") ;

  TFile *bkgFile = TFile::Open(infile.c_str());
  cout <<"p2"<<endl;
  TTree* bkgTree = (TTree*)bkgFile->Get("cattree");
  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,1,1) ;

  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  obsAndWeight.add(*catin);
  obsAndWeight.add(*disc_nw2016);
  obsAndWeight.add(*reco_mu1_eta);
  obsAndWeight.add(*evWeight);

  cout <<"cut: "<<cut<<endl;

  string dataset_ch_name(Form("Sig"+procname+"_cat%d",catindex));
  RooCmdArg cut_arg = cut==""?RooCmdArg::none():RooFit::Cut(cut);
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(),
                RooArgSet(obsAndWeight), RooFit::WeightVar(*evWeight),
                RooFit::Import(*data_chain),
                cut_arg
                );

  //cout << "DEBUG: result of fit on all data " << endl ;
  //RooFitResult* r_full = bkgpdf->fitTo(*data,Save()) ;
  
  
  mdimu->setRange("R1",110,120);
  mdimu->setRange("R2",130,150);

  cout << "DEBUG: result of fit on sideband data " << endl ;
  RooFitResult* r = bkgpdf->fitTo(*data,RooFit::Strategy(1), RooFit::Extended(kTRUE), RooFit::Save(kTRUE), RooFit::Range("R1,R2")) ;
  
  //r_full->Print() ;
  //r->Print() ;

  RooBinning tbins(110,150);
  tbins.addUniform(40,110,150) ;
  RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
  data->plotOn(dtframe,Binning(tbins));
  bkgpdf->plotOn(dtframe);
  TCanvas* c1 = new TCanvas();
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"bkg_13TeV.cat%d_m.png",catindex)));

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;
  w->import(*bkgpdf) ;
  w->import(*data) ;
  w->Print() ;
  w->writeToFile(TString::Format(outfile+"/pdfs/Hmm.input"+procname+"bkg_13TeV.cat%d_m.root",catindex));

  delete data;
  delete bkgpdf;
}

void HiggsMuMuFit::getpdf()
{
  if(HiggsMuMuFit::isbkg.Contains("bkg")) HiggsMuMuFit::bkgfit();
  else HiggsMuMuFit::sigfit();
}
