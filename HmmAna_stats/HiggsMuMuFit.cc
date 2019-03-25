#define HiggsMuMuFit_cxx
#include "HiggsMuMuFit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <RooNLLVar.h>
#include <RooMinimizer.h>
#include <RooNumber.h>
#include <RooGenericPdf.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

//using namespace RooStats;
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
 
  TString proc = isbkgstr;
  //TString sel15 = mass_cut + " cat_index<4"; //ttH
  //TString sel16 = mass_cut + " cat_index>3 && cat_index<7"; //VH-lep
  /*
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
  */
  //if adding VH-lep categories etc mass_cut + " cat_index>6 &&" + 
  bool mass_cut1 =  mass_cut;

  /*2018
  TString sel0 = mass_cut1 + " BDT_incl < -0.427776";

  TString sel1 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > -0.427776 && BDT_incl < 0.0219585";
  TString sel2 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > -0.427776 && BDT_incl < 0.0219585";
  TString sel3 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > -0.427776 && BDT_incl < 0.0219585";

  TString sel4 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.0219585 && BDT_incl < 0.225359";
  TString sel5 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.0219585 && BDT_incl < 0.225359";
  TString sel6 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.0219585 && BDT_incl < 0.225359";

  TString sel7 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.225359 && BDT_incl < 0.372603";
  TString sel8 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.225359 && BDT_incl < 0.372603";
  TString sel9 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.225359 && BDT_incl < 0.372603";

  TString sel10 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.372603 && BDT_incl < 0.613851";
  TString sel11 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.372603 && BDT_incl < 0.613851";
  TString sel12 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.372603 && BDT_incl < 0.613851";

  TString sel13 = mass_cut1 + " BDT_incl > 0.613851 && BDT_incl < 0.717566";
  TString sel14 = mass_cut1 + " BDT_incl > 0.717566";
  */
  /*
  2016
  TString sel0 = mass_cut1 + " BDT_incl < -0.41438";

  TString sel1 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > -0.41438 && BDT_incl < 0.0516541";
  TString sel2 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > -0.41438 && BDT_incl < 0.0516541";
  TString sel3 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > -0.41438 && BDT_incl < 0.0516541";

  TString sel4 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.0516541 && BDT_incl < 0.254417";
  TString sel5 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.0516541 && BDT_incl < 0.254417";
  TString sel6 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.0516541 && BDT_incl < 0.254417";

  TString sel7 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.254417 && BDT_incl < 0.402352";
  TString sel8 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.254417 && BDT_incl < 0.402352";
  TString sel9 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.254417 && BDT_incl < 0.402352";

  TString sel10 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.402352 && BDT_incl < 0.643802";
  TString sel11 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.402352 && BDT_incl < 0.643802";
  TString sel12 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.402352 && BDT_incl < 0.643802";

  TString sel13 = mass_cut1 + " BDT_incl > 0.643802 && BDT_incl < 0.739798";
  TString sel14 = mass_cut1 + " BDT_incl > 0.739798";
 
 2017:
  */
  TString sel0 = mass_cut1 + " BDT_incl < -0.423045";

  TString sel1 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > -0.423045 && BDT_incl < 0.0245217";
  TString sel2 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > -0.423045 && BDT_incl < 0.0245217";
  TString sel3 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > -0.423045 && BDT_incl < 0.0245217";

  TString sel4 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.0245217 && BDT_incl < 0.225491";
  TString sel5 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.0245217 && BDT_incl < 0.225491";
  TString sel6 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.0245217 && BDT_incl < 0.225491";

  TString sel7 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.225491 && BDT_incl < 0.373104";
  TString sel8 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.225491 && BDT_incl < 0.373104";
  TString sel9 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.225491 && BDT_incl < 0.373104";

  TString sel10 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.373104 && BDT_incl < 0.612933";
  TString sel11 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.373104 && BDT_incl < 0.612933";
  TString sel12 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.373104 && BDT_incl < 0.612933";

  TString sel13 = mass_cut1 + " BDT_incl > 0.612933 && BDT_incl < 0.716738";
  TString sel14 = mass_cut1 + " BDT_incl > 0.716738";
  
  /*
  //old 2016 pas?? 
  TString sel0 = mass_cut1 + " BDT_incl < -0.4";

  TString sel1 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > -0.4 && BDT_incl < 0.05";
  TString sel2 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > -0.4 && BDT_incl < 0.05";
  TString sel3 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > -0.4 && BDT_incl < 0.05";

  TString sel4 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.05 && BDT_incl < 0.25";
  TString sel5 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.05 && BDT_incl < 0.25";
  TString sel6 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.05 && BDT_incl < 0.25";

  TString sel7 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.25 && BDT_incl < 0.40";
  TString sel8 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.25 && BDT_incl < 0.40";
  TString sel9 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.25 && BDT_incl < 0.40";

  TString sel10 = mass_cut1 + " fabs(max_reco_mu_eta)<0.9 && BDT_incl > 0.40 && BDT_incl < 0.65";
  TString sel11 = mass_cut1 + " fabs(max_reco_mu_eta)>0.9 && fabs(max_reco_mu_eta)<1.9  && BDT_incl > 0.40 && BDT_incl < 0.65";
  TString sel12 = mass_cut1 + " fabs(max_reco_mu_eta)>1.9 && BDT_incl > 0.40 && BDT_incl < 0.65";

  TString sel13 = mass_cut1 + " BDT_incl > 0.65 && BDT_incl < 0.73";
  TString sel14 = mass_cut1 + " BDT_incl > 0.73";
  */ 
  //vector<TString> selection = {sel0, sel1, sel2, sel3, sel4, sel5, sel6, sel7, sel8, sel9, sel10, sel11, sel12, sel13, sel14};

  vector<TString> selection = {sel5};

  int category = selection.size();
  HiggsMuMuFit hmm(inputFileList, outFileName, isbkgstr, selection, proc, category);

  //hmm.set(inputFileList, outFileName, isbkgstr, selection, proc, category);
  hmm.getpdf();

  return 0;
}

void HiggsMuMuFit::sigfit_DCB(){

  /*
  cout <<"sigfit: double-sided CB"<<endl;

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  RooRealVar* catin= new RooRealVar("cat_index","cat_index",1,0,100) ;
  RooRealVar* max_reco_mu_eta = new RooRealVar("max_reco_mu_eta","max_reco_mu_eta",1,-100,100) ;
  RooRealVar* BDT_incl = new RooRealVar("BDT_incl","BDT_incl",0,-1,1) ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",110,110,150) ;
  mdimu->setBins(800);
  mdimu->setRange("RF",110,150);
  RooRealVar* cut_based_ct = new RooRealVar("cat_index","cat_index",1,0,10) ;

 for(int catindex=0; catindex<ncat; catindex++){
  RooRealVar* mean = new RooRealVar(Form(procname+"cat%d_mean",catindex),Form("mean for "+procname+"cat%d", catindex),125, 110, 145) ;
  RooRealVar* sigma = new RooRealVar(Form(procname+"cat%d_sigma",catindex),Form("sigma for "+procname+"cat%d", catindex),1, 0.1, 50) ;
  RooRealVar* alpha1 = new RooRealVar(Form(procname+"cat%d_alpha1",catindex),Form("alpha1 for "+procname+"cat%d", catindex),1, 0.05, 50) ;
  RooRealVar* alpha2 = new RooRealVar(Form(procname+"cat%d_alpha2",catindex),Form("alpha2 for "+procname+"cat%d", catindex),1, 0.05, 50) ;
  RooRealVar* n1 = new RooRealVar(Form(procname+"cat%d_n1",catindex),Form("n1 for "+procname+"cat%d", catindex),2, 0.05, 50) ;
  RooRealVar* n2 = new RooRealVar(Form(procname+"cat%d_n2",catindex),Form("n2 for "+procname+"cat%d", catindex),2, 0.05, 50) ;
  RooDoubleCB* dimupdf = new RooDoubleCB(Form(procname+"cat%d_sig_pdf",catindex),Form(procname+"cat%d_sig_pdf",catindex), *mdimu, *mean, *sigma, *alpha1, *n1, *alpha2, *n2);

  RooRealVar* nsig = new RooRealVar(Form(procname+"cat%d_nsig",catindex),Form("number of signal events in signalRange "+procname+"cat%d", catindex),10000,0,100000000) ;
  RooExtendPdf* sigpdf = new RooExtendPdf(Form(procname+"cat%d_sig_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*dimupdf,*nsig,"RF") ;

  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-100000,100000) ;

  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  obsAndWeight.add(*catin);
  obsAndWeight.add(*BDT_incl);
  obsAndWeight.add(*max_reco_mu_eta);
  obsAndWeight.add(*evWeight);

  TString cut = cuts[catindex];
  cout <<"cut: "<<cut<<endl;
  string dataset_ch_name(Form("Sig"+procname+"_cat%d",catindex));
  RooCmdArg cut_arg = cut==""?RooCmdArg::none():RooFit::Cut(cut);
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(),
                RooArgSet(obsAndWeight), RooFit::WeightVar(*evWeight),
                RooFit::Import(*data_chain),
                cut_arg
                );

  RooNLLVar* nll = (RooNLLVar*)sigpdf->createNLL(*data, Extended(kTRUE));
  nll->enableOffsetting(true);
  RooMinimizer minim(*nll);
  minim.setStrategy(1);
  minim.setPrintLevel(0);
  int status = minim.minimize("Minuit2", "Migrad");
  if(status!=0){
      cout <<"fit to data does not converge, try again"<<endl;
      RooMinimizer minim2(*nll);
      minim2.setStrategy(1);
      minim2.setPrintLevel(0);
      int status2 = minim2.minimize("Minuit2", "Migrad");
      if(status2!=0){ 
        cout <<"fit to data does not converge, try again"<<endl;
        RooMinimizer minim3(*nll);
        minim3.setStrategy(1);
        minim3.setPrintLevel(0);
        int status3 = minim3.minimize("Minuit2", "Migrad");
        if(status3!=0){
            cout <<"fit to data does not converge, after three fits"<<endl;
        } 
      } 
  }   
  
  RooBinning tbins(110,150);
  tbins.addUniform(80,110,150) ;
  RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
  data->plotOn(dtframe,Binning(tbins));
  sigpdf->plotOn(dtframe);

  sigpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(0))) ;
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextSize(0.03);

  TCanvas* c1 = new TCanvas();
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m.png",catindex)));

  gPad->SetLogy();
  double maxv = dtframe->GetMaximum();
  dtframe->SetMinimum(0.01*maxv);
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m_log.png",catindex)));

  nsig->setConstant(true);
  mean->setConstant(true);
  sigma->setConstant(true);
  n1->setConstant(true);
  n2->setConstant(true);
  alpha1->setConstant(true);
  alpha2->setConstant(true);

  w->import(*sigpdf) ;
  w->import(*data) ;

  delete data_chain;
  delete data;
  delete sigpdf;
  }
  w->Print() ;
  w->writeToFile(outfile+"/pdfs/Hmm.input"+procname+"sig_13TeV.root");
  */
}

void HiggsMuMuFit::sigfit()
{ 

  cout <<"sigfit: sum of up to three Gaussian functions"<<endl;

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  RooRealVar* catin= new RooRealVar("cat_index","cat_index",1,0,100) ;
  RooRealVar* max_reco_mu_eta = new RooRealVar("max_reco_mu_eta","max_reco_mu_eta",1,-100,100) ;
  RooRealVar* BDT_incl = new RooRealVar("BDT_incl","BDT_incl",0,-1,1) ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",110,110,150) ;
  mdimu->setBins(800);
  mdimu->setRange("RF",110,150);
  RooRealVar* cut_based_ct = new RooRealVar("cat_index","cat_index",1,0,10) ;

 for(int catindex=0; catindex<ncat; catindex++){
  RooRealVar* mean1 = new RooRealVar(Form(procname+"cat%d_mean1",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
  RooRealVar* mean2 = new RooRealVar(Form(procname+"cat%d_mean2",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
  RooRealVar* mean3 = new RooRealVar(Form(procname+"cat%d_mean3",catindex),Form("mean3 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;

  RooRealVar* sigma1 = new RooRealVar(Form(procname+"cat%d_sigma1",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex),3, 0.0, 10) ;
  RooRealVar* sigma2 = new RooRealVar(Form(procname+"cat%d_sigma2",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex),5, 0.0, 40) ;
  RooRealVar* sigma3 = new RooRealVar(Form(procname+"cat%d_sigma3",catindex),Form("sigma3 of gaussians for "+procname+"cat%d", catindex),10, 0.0, 40) ;

  RooGaussian* sig1 = new RooGaussian(Form(procname+"cat%d_sig1",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
  RooGaussian* sig2 = new RooGaussian(Form(procname+"cat%d_sig2",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;
  RooGaussian* sig3 = new RooGaussian(Form(procname+"cat%d_sig3",catindex),Form("Signal component 3 for "+procname+"cat%d", catindex),*mdimu, *mean3, *sigma3) ;

  RooRealVar* sig1frac = new RooRealVar(Form(procname+"cat%d_sig1frac",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),0.8,0.,1.) ;
  RooRealVar* sig2frac = new RooRealVar(Form(procname+"cat%d_sig2frac",catindex),Form("signal fraction of component 2 for "+procname+"cat%d", catindex),0.1,0.,1.) ;
  RooAddPdf* dimupdf = new RooAddPdf(Form(procname+"cat%d_sig_pdf",catindex),Form(procname+"cat%d_sig_pdf",catindex),RooArgList(*sig1,*sig2,*sig3),RooArgList(*sig1frac,*sig2frac)) ;
  //RooAddPdf* dimupdf = new RooAddPdf(Form(procname+"cat%d_sig_pdf",catindex),Form(procname+"cat%d_sig_pdf",catindex),RooArgList(*sig1,*sig2),RooArgList(*sig1frac)) ;

  RooRealVar* nsig = new RooRealVar(Form(procname+"cat%d_nsig",catindex),Form("number of signal events in signalRange from fit "+procname+"cat%d", catindex),10000,0,100000000) ;
  RooExtendPdf* sigpdf = new RooExtendPdf(Form(procname+"cat%d_sig_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*dimupdf,*nsig,"RF") ;

  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-100000,100000) ;

  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  obsAndWeight.add(*catin);
  obsAndWeight.add(*BDT_incl);
  obsAndWeight.add(*max_reco_mu_eta);
  obsAndWeight.add(*evWeight);

  TString cut = cuts[catindex];
  cout <<"cut: "<<cut<<endl;
  string dataset_ch_name(Form("Sig"+procname+"_cat%d",catindex));
  RooCmdArg cut_arg = cut==""?RooCmdArg::none():RooFit::Cut(cut);
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(),
                RooArgSet(obsAndWeight), RooFit::WeightVar(*evWeight),
                RooFit::Import(*data_chain),
                cut_arg
                );

  double sum = data->sumEntries();
  if(procname.Contains("ggH")){
     sum = sum*1.100589035;
     cout <<"scale to N3LO xs"<<endl;
  }
  RooRealVar* nsig_sum = new RooRealVar(Form(procname+"cat%d_nsig_sum",catindex),Form("number of signal events in signalRange "+procname+"cat%d", catindex),sum) ;
  nsig_sum->setConstant(true);

  //RooFitResult* r = sigpdf->fitTo(*data,RooFit::Strategy(1), RooFit::Extended(kTRUE), RooFit::Save(kTRUE));
  //r->Print() ;

  RooNLLVar* nll = (RooNLLVar*)sigpdf->createNLL(*data, Extended(kTRUE));
  nll->enableOffsetting(true);
  RooMinimizer minim(*nll);
  minim.setStrategy(1);
  minim.setPrintLevel(0);
  int status = minim.minimize("Minuit2", "Migrad");
  if(status!=0){
      cout <<"fit to data does not converge, try again"<<endl;
      RooMinimizer minim2(*nll);
      minim2.setStrategy(1);
      minim2.setPrintLevel(0);
      int status2 = minim2.minimize("Minuit2", "Migrad"); 
      if(status2!=0){ 
        cout <<"fit to data does not converge, try again"<<endl;
        RooMinimizer minim3(*nll);
        minim3.setStrategy(1);
        minim3.setPrintLevel(0);
        int status3 = minim3.minimize("Minuit2", "Migrad");
        /*
        if(status3!=0){ 
          cout <<"fit to data does not converge, reduce to two Gaussians, try again"<<endl;
          sig2frac->setVal(0);
          sig2frac->setConstant(true);
          sigma2->setConstant(true);
          mean2->setConstant(true);
          RooMinimizer minim4(*nll);
          minim4.setStrategy(1);
          minim4.setPrintLevel(0);
          int status4 = minim4.minimize("Minuit2", "Migrad");
          if(status4!=0){
            cout <<"fit to data does not converge, after four fits"<<endl;
          } 
        }
        */
      }
  }

  cout <<"nsig_sum: "<<sum<<endl;
  nsig->Print();

  RooBinning tbins(110,150);
  tbins.addUniform(80,110,150) ;
  RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
  data->plotOn(dtframe,Binning(tbins));
  sigpdf->plotOn(dtframe);
  sigpdf->plotOn(dtframe, Components(*sig1), LineColor(kRed), LineStyle(kDashed));
  sigpdf->plotOn(dtframe, Components(*sig2), LineColor(kGreen), LineStyle(kDashed));
  sigpdf->plotOn(dtframe, Components(*sig3), LineColor(kOrange), LineStyle(kDashed));

  //data->statOn(dtframe,Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(0))) ;
  sigpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(0))) ;
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextSize(0.03);

  TCanvas* c1 = new TCanvas();
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m.png",catindex)));

  gPad->SetLogy();
  double maxv = dtframe->GetMaximum();
  dtframe->SetMinimum(0.01*maxv);
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m_log.png",catindex)));

  if((sig1frac->getVal()+sig2frac->getVal())>1.0){
    if(sig1frac->getVal() > sig2frac->getVal()) sig2frac->setVal(0.);
    else sig1frac->setVal(0.);
  }
  nsig->setConstant(true);
  sig1frac->setConstant(true);
  sig2frac->setConstant(true);
  sigma1->setConstant(true);
  sigma2->setConstant(true);
  sigma3->setConstant(true);
  mean1->setConstant(true);
  mean2->setConstant(true);
  mean3->setConstant(true);

  w->import(*nsig_sum);
  w->import(*sigpdf);
  w->import(*data);

  delete data_chain;
  delete data;
  delete sigpdf; 
  }
  w->Print() ;
  w->writeToFile(outfile+"/pdfs/Hmm.input"+procname+"sig_13TeV.root");
}   

void HiggsMuMuFit::bkgfit()
{
}

void HiggsMuMuFit::bkgfit_BWZReduxMBern()
{
  cout <<"bkgfit: BWZRedu x MBern"<<endl;

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",125,110,150) ;
  mdimu->setRange("R1",110,120);
  mdimu->setRange("R2",130,150);
  mdimu->setRange("RF",110,150);
  mdimu->setBins(800);
  RooFormulaVar* mdimu_range =  new RooFormulaVar("mdimu_range","(@0-110.0)/(40.0)",RooArgSet(*mdimu));

 for(int catindex=0; catindex<ncat; catindex++){
  RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
  RooRealVar* btwzr_a2 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_btwzr_a2_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_btwzr_a2_cat%d",catindex),0.1,-20,20) ;
  RooRealVar* btwzr_a3 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_btwzr_a3_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_btwzr_a3_cat%d",catindex),0.1,-20,20) ;
  RooRealVar* b1 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;
  RooRealVar* b2 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_b2_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_b2_cat%d",catindex),0.01,-1.0,1.0) ;
  RooRealVar* b3 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_b3_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_b3_cat%d",catindex),0.01,-1.0,1.0) ;
  RooRealVar* b4 = new RooRealVar(TString::Format("CMS_dimu_"+procname+"_b4_cat%d",catindex),TString::Format("CMS_dimu_"+procname+"_b4_cat%d",catindex),0.01,-1.0,1.0) ;
  RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3));

  //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*5*pow(@4,4)+20*@6*pow(@4,3)*(1-@4)+@7*30*@4*@4*(1-@4)*(1-@4) + @8 *20*@4*pow(1-@4,3)+(-@5-@6-@7-@8)*5*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4));
  //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*@5*5*pow(@4,4)+20*@6*@6*pow(@4,3)*(1-@4)+@7*@7*30*@4*@4*(1-@4)*(1-@4) + @8*@8*20*@4*pow(1-@4,3)+(-@5*5-@6*6-@7*@7-@8*@8)*5*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4)); 
 //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * ((1-@5*@5-@6*@6-@7*@7-@8*@8)*pow(@4,4)+@5*@5*4*pow(@4,3)*(1-@4)+@6*@6*6*@4*@4*(1-@4)*(1-@4)+@7*@7*4*@4*pow(1-@4,3)+@8*@8*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4));

  RooRealVar* nbkg = new RooRealVar(Form(procname+"cat%d_nbkg",catindex),Form("number of background events in signalRange "+procname+"cat%d", catindex),500,-RooNumber::infinity(),RooNumber::infinity()) ;
  RooExtendPdf* bkgpdf = new RooExtendPdf(Form(procname+"cat%d_bkg_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*mdimupdf,*nbkg,"RF") ;

  //TFile *bkgFile = TFile::Open(infile.c_str());
  cout <<"p2"<<endl;
  //TTree* bkgTree = (TTree*)bkgFile->Get("cattree");
  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,1,1) ;

  bool weighted = false;
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  if(weighted) obsAndWeight.add(*evWeight);

  TString cut = cuts[catindex];
  cout <<"cut: "<<cut<<endl;
  RooCmdArg wgt_arg = weighted ? RooFit::WeightVar("evt_weight") : RooCmdArg::none() ;

  TTree* cutChain = data_chain->CopyTree(cut);
  TTree* cutChainblind = data_chain->CopyTree(cut+" && (Higgs_mass>130 || Higgs_mass<120)");
  string dataset_ch_name(Form("Sig"+procname+"_cat%d",catindex));
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain), wgt_arg);
  RooDataSet* datablind = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChainblind), wgt_arg);
  cout << "DEBUG: result of fit on sideband data " << endl ;
  RooNLLVar* nll_fix = (RooNLLVar*)bkgpdf->createNLL(*data, Extended(kTRUE), Range("RF"));
  RooMinimizer minim_fix(*nll_fix);
  minim_fix.setStrategy(1);
  minim_fix.setPrintLevel(0);
  int status_fix = minim_fix.minimize("Minuit2", "Migrad");
  if(status_fix!=0) cout <<"fit to data does not converge"<<endl;
  //btwzr_a1->setConstant(true);
  //btwzr_a2->setConstant(true);
  //btwzr_a3->setConstant(true);
  cout <<"1. expected events: "<<bkgpdf->expectedEvents(*mdimu)<<endl;
  nbkg->setVal(bkgpdf->expectedEvents(*mdimu));

  RooBinning tbins(110,150);
  tbins.addUniform(80,110,150) ;
  RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
  datablind->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(0)));
  bkgpdf->plotOn(dtframe,NormRange("R1,R2"),Range("RF"));

  bkgpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(0))) ;
  dtframe->getAttText(Form(procname+"cat%d_bkg_extpdf_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_bkg_extpdf_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_bkg_extpdf_paramBox",catindex))->SetTextSize(0.03);
  cout <<"2. expected events: "<<bkgpdf->expectedEvents(*mdimu)<<endl;

  TCanvas* c1 = new TCanvas();
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"bkg_13TeV.cat%d_m.png",catindex)));

  gPad->SetLogy();
  double maxv = dtframe->GetMaximum();
  dtframe->SetMinimum(0.01*maxv);
  dtframe->Draw();
  c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"bkg_13TeV.cat%d_m_log.png",catindex)));

  w->import(*bkgpdf, RecycleConflictNodes()) ;
  w->import(*data) ;
  delete data;
  delete bkgpdf;

 }
  w->Print() ;
  w->writeToFile(outfile+"/pdfs/Hmm.input"+procname+"bkg_13TeV.root");
}

void HiggsMuMuFit::getpdf()
{
  if(HiggsMuMuFit::isbkg.Contains("bkg")) HiggsMuMuFit::bkgfit_BWZReduxMBern(); //bkgfit();
  else HiggsMuMuFit::sigfit();
}
