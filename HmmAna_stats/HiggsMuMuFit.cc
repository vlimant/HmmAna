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
  TString funcname = argv[4];
  TString useweight = argv[5];

  TString mass_cut = "Higgs_mass > 110 && Higgs_mass< 150 &&";
 
  TString proc = isbkgstr;
  TString mass_cut1 =  mass_cut;

  //TString sel0 = mass_cut1 + " disc_simpleNN > 0.1"; //disc_advNN
  //TString sel1 = mass_cut1 + " disc_simpleNN < 0.1 && disc_simpleNN>0.0"; 

  TString sel0 = mass_cut1 + " disc_advNN > 0.1";
  TString sel1 = mass_cut1 + " disc_advNN < 0.1 && disc_advNN>0.0"; 

  vector<TString> selection = {sel0,sel1};

  int category = selection.size();
  
  HiggsMuMuFit hmm(inputFileList, outFileName, isbkgstr, funcname, selection, proc, category, useweight);
  Utils::setDefaultMinimize();
  hmm.getpdf();

  return 0;
}

void HiggsMuMuFit::sigfit_DCB(){
   //to be added
}

void HiggsMuMuFit::sigfit()
{ 

  cout <<"sigfit: sum of up to three Gaussian functions"<<endl;

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  //RooRealVar* catin= new RooRealVar("cat_index","cat_index",1,0,100) ;
  //RooRealVar* max_reco_mu_eta = new RooRealVar("max_reco_mu_eta","max_reco_mu_eta",1,-100,100) ;
  //RooRealVar* BDT_incl = new RooRealVar("BDT_incl","BDT_incl",0,-1,1) ;
  
  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",110,110,150) ;
  mdimu->setBins(800);
  mdimu->setRange("RF",110,150);
  //RooRealVar* cut_based_ct = new RooRealVar("cat_index","cat_index",1,0,10) ;

 for(int catindex=0; catindex<ncat; catindex++){
  
  RooRealVar* mean1 = new RooRealVar(Form("mean1_"+procname+"cat%d",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
  RooRealVar* mean2 = new RooRealVar(Form("mean2_"+procname+"cat%d",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
  RooRealVar* mean3 = new RooRealVar(Form("mean3_"+procname+"cat%d",catindex),Form("mean3 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;

  RooRealVar* sigma1 = new RooRealVar(Form("sigma1_"+procname+"cat%d",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex),3, 0.0, 10) ;
  RooRealVar* sigma2 = new RooRealVar(Form("sigma2_"+procname+"cat%d",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex),5, 0.0, 40) ;
  RooRealVar* sigma3 = new RooRealVar(Form("sigma3_"+procname+"cat%d",catindex),Form("sigma3 of gaussians for "+procname+"cat%d", catindex),10, 0.0, 40) ;

  RooGaussian* sig1 = new RooGaussian(Form("sig1_"+procname+"cat%d",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
  RooGaussian* sig2 = new RooGaussian(Form("sig2_"+procname+"cat%d",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;
  RooGaussian* sig3 = new RooGaussian(Form("sig3_"+procname+"cat%d",catindex),Form("Signal component 3 for "+procname+"cat%d", catindex),*mdimu, *mean3, *sigma3) ;

  RooRealVar* sig1frac = new RooRealVar(Form("sig1frac_"+procname+"cat%d",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),0.8,0.,1.) ;
  RooRealVar* sig2frac = new RooRealVar(Form("sig1frac_"+procname+"cat%d",catindex),Form("signal fraction of component 2 for "+procname+"cat%d", catindex),0.1,0.,1.) ;
  RooAddPdf* dimupdf = new RooAddPdf(Form("sig_pdf_"+procname+"cat%d",catindex),Form(procname+"cat%d_sig_pdf",catindex),RooArgList(*sig1,*sig2),RooArgList(*sig1frac)) ;
  
  //RooAbsPdf* dimupdf = HiggsMuMuFit::SGauss( mdimu, catindex);
  RooRealVar* nsig = new RooRealVar(Form(procname+"cat%d_nsig",catindex),Form("number of signal events in signalRange from fit "+procname+"cat%d", catindex),10000,0,100000000) ;
  RooExtendPdf* sigpdf = new RooExtendPdf(Form(procname+"cat%d_sig_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*dimupdf,*nsig,"RF") ;

  TChain* data_chain = loader(infile.c_str(), "cattree");

  RooRealVar* evWeight = new RooRealVar("genWeight","genWeight",1,-100000,100000) ;

  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  //obsAndWeight.add(*catin);
  //obsAndWeight.add(*BDT_incl);
  //obsAndWeight.add(*max_reco_mu_eta);
  obsAndWeight.add(*evWeight);

  TString cut = cuts[catindex];
  cout <<"cut: "<<cut<<endl;
  TTree* cutChain = data_chain->CopyTree(cut);
  string dataset_ch_name(Form("Sig_"+procname+"_cat%d",catindex));
  RooCmdArg wgt_arg = RooFit::WeightVar("genWeight");
  RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain), wgt_arg);
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
  //sigpdf->plotOn(dtframe, Components(*sig3), LineColor(kOrange), LineStyle(kDashed));

  //data->statOn(dtframe,Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(0))) ;
  sigpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(2))) ;
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form(procname+"cat%d_sig_extpdf_paramBox",catindex))->SetTextSize(0.03);

  TCanvas* c1 = new TCanvas("c1","c1");
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

void HiggsMuMuFit::sbfit()
{
  if(funclist.size()!=ncat){
    cout <<"error "<<funclist.size()<<"vs"<<ncat<<endl;
    return;
  }

  cout <<"signal + bkg fit in functions: "<<endl;
  for(int i=0; i<funclist.size(); i++) cout <<" cat"<<i<<":"<<funclist.at(i);
  cout <<endl;

  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",125,110,150) ;
  mdimu->setRange("R1",110,120);
  mdimu->setRange("R2",130,150);
  mdimu->setRange("RF",110,150);
  mdimu->setBins(800);

  RooFormulaVar* mdimu_range =  new RooFormulaVar("mdimu_range","(@0-110.0)/(40.0)",RooArgSet(*mdimu));
  RooRealVar* evWeight = new RooRealVar("genWeight","genWeight",1,-100000,100000) ;

  TChain* data_chain = loader(infile.c_str(), "cattree");
  TChain* data_chain_sig = loader("/eos/cms/store/user/nlu/Hmm/VBF/categorization/VBFHToMuMu_2018_NNscore.root", "cattree");
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  if(usew.EqualTo("T")) obsAndWeight.add(*evWeight);

  for(int catindex=0; catindex<ncat; catindex++){
        TString cut = cuts[catindex];
        cout <<"cut: "<<cut<<endl;
        RooCmdArg wgt_arg = usew ? RooFit::WeightVar("genWeight") : RooCmdArg::none() ;

        TTree* cutChain_sig = data_chain_sig->CopyTree(cut);
        string dataset_ch_name_sig(Form("sig_"+procname+"_cat%d",catindex));
        RooDataSet* datasig = new RooDataSet(dataset_ch_name_sig.c_str(), dataset_ch_name_sig.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain_sig), wgt_arg);

        RooAbsPdf* bkgpdf = HiggsMuMuFit::bkgfunc( mdimu, mdimu_range, catindex, funclist.at(catindex));
        RooAbsPdf* sigpdf = HiggsMuMuFit::SGauss( mdimu, catindex, 124.845, 121.643, 1.7824, 5.80949, 0.819765);
        RooRealVar* nsigexp = new RooRealVar(Form("nsigexp_"+procname+"cat%d",catindex),"no. expected signal events",datasig->sumEntries(),datasig->sumEntries(),datasig->sumEntries());
        nsigexp->setConstant(true);
        RooRealVar* mu = new RooRealVar(Form("mu_"+procname+"cat%d",catindex),"signal strength",1.,0.0,1000.0);
        RooFormulaVar* nsig =  new RooFormulaVar(Form("nsig_"+procname+"cat%d",catindex),"@0*@1",RooArgSet(*nsigexp,*mu));
        RooExtendPdf* smodel = new RooExtendPdf(Form("smodel_"+procname+"cat%d",catindex),"extended signal p.d.f",*sigpdf,*nsig,"RF");
        RooRealVar* nbkg = new RooRealVar(Form("nbkg_"+procname+"cat%d",catindex),"no. bkg ",10.,0.0,RooNumber::infinity());
        RooExtendPdf* bmodel = new RooExtendPdf(Form("bmodel_"+procname+"cat%d",catindex),"extended bkg p.d.f",*bkgpdf,*nbkg,"RF");
        //RooRealVar* sigfrac = new RooRealVar(Form("sigfrac_"+procname+"cat%d",catindex),"sig frac",0.5);
        //sigfrac->setConstant(true);
        RooAbsPdf* model = new RooAddPdf(Form("sbmodel_"+procname+"cat%d",catindex),"signal + bkg",RooArgList(*smodel,*bmodel)) ;

        TTree* cutChain = data_chain->CopyTree(cut);
        TTree* cutChainblind = data_chain->CopyTree(cut+" && (Higgs_mass>130 || Higgs_mass<120)");
        string dataset_ch_name(Form("data_"+procname+"_cat%d",catindex));
        RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain), wgt_arg);
        RooDataSet* datablind = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChainblind), wgt_arg);

        cout << "DEBUG: result of fit on full data " << endl ;
        RooNLLVar* nll_fix = (RooNLLVar*)model->createNLL(*data, Extended(kTRUE), Range("RF"));
        int status_fix = Utils::minimizeMinosTest(nll_fix,mu,1.0); //Utils::minimizeTest(nll_fix,1.0);
        if(status_fix!=0){
            cout <<"fit to data does not converge"<<endl;
            status_fix = Utils::minimizeMinosTest(nll_fix,mu,1.0);
            if(status_fix!=0){ 
                   cout <<"fit to data does not converge after a second trial"<<endl;
            }
        }

        InitTreeVars();
        t_category = catindex;
        t_mu = mu->getVal();
        t_errlo = mu->getErrorLo();
        t_errhi = mu->getErrorHi();
        t_err = 0.5*(mu->getErrorHi()-mu->getErrorLo());
        t_nsig = smodel->expectedEvents(*mdimu);
        t_nbkg = bmodel->expectedEvents(*mdimu);
        outtree->Fill();

        cout <<"1. expected events from fit: "<<model->expectedEvents(*mdimu)<<endl;

        RooBinning tbins(110,150);
        tbins.addUniform(80,110,150) ;
        RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
        //dtframe->getAttText()->SetTextSize(0.025) ;
        //if(!usew.EqualTo("T")) datablind->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
        //else{
        data->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
        cout <<"expected number of events from MC: " <<data->sumEntries()<<endl;
        //}
        //model->plotOn(dtframe,NormRange("R1,R2"),Range("RF"));
        model->plotOn(dtframe, Normalization(1.0,RooAbsReal::RelativeExpected), Range("RF"));
        //model->plotOn(dtframe, Components(*smodel), Normalization(1.0,RooAbsReal::RelativeExpected), LineColor(kRed), LineStyle(kDashed));
        //model->plotOn(dtframe, Components(*bmodel), Normalization(1.0,RooAbsReal::RelativeExpected), LineColor(kGreen), LineStyle(kDashed));
        model->paramOn(dtframe,Layout(0.25, 0.75, 0.9), Format("NEU",AutoPrecision(3))) ;
        smodel->plotOn(dtframe, Normalization(1.0,RooAbsReal::RelativeExpected), Range("RF"),LineColor(kRed), LineStyle(kDashed));
        cout <<"err2"<<endl;
        bmodel->plotOn(dtframe, Normalization(1.0,RooAbsReal::RelativeExpected), Range("RF"),LineColor(kGreen), LineStyle(kDashed));
        cout <<"err3"<<endl;
        //dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextColor(kBlack);
        //dtframe->getAttFill(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetFillStyle(0);
        //dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextSize(0.03);
        cout <<"2. expected events from fit: "<<model->expectedEvents(*mdimu)<<endl;
        cout <<"2. expected signal events from fit: "<<smodel->expectedEvents(*mdimu)<<" vs MC "<<datasig->sumEntries()<<endl;
        cout <<"2. expected bkg events from fit: "<<bmodel->expectedEvents(*mdimu)<<" vs MC "<<data->sumEntries()-datasig->sumEntries()<<endl;
        TCanvas* c1 = new TCanvas("c1","c1");
        double maxv = dtframe->GetMaximum();
        dtframe->SetMaximum(1.5*maxv);
        dtframe->Draw();
        c1->SaveAs(Form(outfile+"/plots/Hmm."+procname+"sigbkg_13TeV.cat%d_m.png",catindex));

        gPad->SetLogy();
        dtframe->SetMaximum(200.*maxv);
        dtframe->SetMinimum(0.01*maxv);
        dtframe->Draw();
        c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sigbkg_13TeV.cat%d_m_log.png",catindex)));
        RooRealVar* Nevt_mc = new RooRealVar("Nevt_mc","Nevt_mc",datasig->sumEntries()) ;
        //RooRealVar* Nevt_sigfit = new RooRealVar("Nevt_sigfit","Nevt_sigfit",nevents->getVal()*(1.0-bkgfrac->getVal())) ;
        //RooRealVar* Nevt_bkgfit = new RooRealVar("Nevt_bkgfit","Nevt_bkgfit",nevents->getVal()*bkgfrac->getVal()) ;

        w->import(*model, RecycleConflictNodes()) ;
        w->import(*data);
        w->import(*datasig);
        w->import(*Nevt_mc);
        //w->import(*Nevt_sigfit);
        //w->import(*Nevt_bkgfit);
        delete data;
        delete model;
        delete cutChain_sig;

    }
    w->Print() ;
    w->writeToFile(outfile+"/pdfs/Hmm.input"+procname+"sigbkg_13TeV.root");
    delete data_chain;
}

void HiggsMuMuFit::bkgfit()
{
  if(funclist.size()!=ncat){
    cout <<"error "<<funclist.size()<<"vs"<<ncat<<endl;
    return;
  }

  cout <<"bkgfit in functions: "<<endl; 
  for(int i=0; i<funclist.size(); i++) cout <<" cat"<<i<<":"<<funclist.at(i);
  cout <<endl;
 
  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;

  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",125,110,150) ;
  mdimu->setRange("R1",110,120);
  mdimu->setRange("R2",130,150);
  mdimu->setRange("RF",110,150);
  mdimu->setBins(800);
  RooFormulaVar* mdimu_range =  new RooFormulaVar("mdimu_range","(@0-110.0)/(40.0)",RooArgSet(*mdimu));
  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-100000,100000) ;
   
  TChain* data_chain = loader(infile.c_str(), "cattree");
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  if(usew.EqualTo("T")) obsAndWeight.add(*evWeight);

  for(int catindex=0; catindex<ncat; catindex++){

        RooAbsPdf* mdimupdf = NULL;
        if(funclist.at(catindex).EqualTo("BWZRedux")){
            mdimupdf = HiggsMuMuFit::BWZRedux( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZReduxB4")){
            mdimupdf = HiggsMuMuFit::BWZReduxB4( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZReduxB3")){
            mdimupdf = HiggsMuMuFit::BWZReduxB3( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZ")){
            mdimupdf = HiggsMuMuFit::BWZ( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZPol1")){
            mdimupdf = HiggsMuMuFit::BWZPol1( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZPol2")){
            mdimupdf = HiggsMuMuFit::BWZPol2( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("BWZGamma")){
            mdimupdf = HiggsMuMuFit::BWZGamma( mdimu, mdimu_range, catindex);
        }
        else if(funclist.at(catindex).EqualTo("SExp")){
            mdimupdf = HiggsMuMuFit::SExp( mdimu, mdimu_range, catindex);
        }
        //else if(funclist.at(catindex).EqualTo("All")){
            //bkgpdf = HiggsMuMuFit::SExp( mdimu, mdimu_range, catindex);
        //}
        else{
          cout <<"function does not find, error"<<endl;
          continue;
        }
     
        RooRealVar* nbkg = new RooRealVar(Form("nbkg_"+procname+"cat%d",catindex),Form("number of background events in signalRange "+procname+"cat%d", catindex),500,0.0,RooNumber::infinity()) ;
        RooExtendPdf* bkgpdf = new RooExtendPdf(Form("extpdf_"+procname+"cat%d_bkg",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*mdimupdf,*nbkg,"RF") ;

        TString cut = cuts[catindex];
        cout <<"cut: "<<cut<<endl;
        RooCmdArg wgt_arg = usew ? RooFit::WeightVar("evt_weight") : RooCmdArg::none() ;

        TTree* cutChain = data_chain->CopyTree(cut);
        TTree* cutChainblind = data_chain->CopyTree(cut+" && (Higgs_mass>130 || Higgs_mass<120)");
        string dataset_ch_name(Form("data_"+procname+"_cat%d",catindex));
        RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain), wgt_arg);
        RooDataSet* datablind = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChainblind), wgt_arg);

        cout << "DEBUG: result of fit on full data " << endl ;
        RooNLLVar* nll_fix = (RooNLLVar*)bkgpdf->createNLL(*data, Extended(kTRUE), Range("RF"));
        RooMinimizer minim_fix(*nll_fix);
        minim_fix.setStrategy(1);
        minim_fix.setPrintLevel(0);
        int status_fix = minim_fix.minimize("Minuit2", "Migrad");
        if(status_fix!=0) cout <<"fit to data does not converge"<<endl;

        cout <<"1. expected events: "<<bkgpdf->expectedEvents(*mdimu)<<endl;
        //nbkg->setVal(bkgpdf->expectedEvents(*mdimu));

        RooBinning tbins(110,150);
        tbins.addUniform(80,110,150) ;
        RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
        if(!usew.EqualTo("T")) datablind->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
        else{
           data->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
           cout <<"expected number of events: " <<data->sumEntries()<<endl;
        }
        bkgpdf->plotOn(dtframe,NormRange("R1,R2"),Range("RF"));

        bkgpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(2))) ;
        dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextColor(kBlack);
        dtframe->getAttFill(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetFillStyle(0);
        dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextSize(0.03);
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
    delete data_chain;
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
  RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3));

  //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*5*pow(@4,4)+20*@6*pow(@4,3)*(1-@4)+@7*30*@4*@4*(1-@4)*(1-@4) + @8 *20*@4*pow(1-@4,3)+(-@5-@6-@7-@8)*5*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4));
  //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*@5*5*pow(@4,4)+20*@6*@6*pow(@4,3)*(1-@4)+@7*@7*30*@4*@4*(1-@4)*(1-@4) + @8*@8*20*@4*pow(1-@4,3)+(-@5*5-@6*6-@7*@7-@8*@8)*5*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4)); 
 //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format(procname+"cat%d_bkg_pdf",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * ((1-@5*@5-@6*@6-@7*@7-@8*@8)*pow(@4,4)+@5*@5*4*pow(@4,3)*(1-@4)+@6*@6*6*@4*@4*(1-@4)*(1-@4)+@7*@7*4*@4*pow(1-@4,3)+@8*@8*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4));

  RooRealVar* nbkg = new RooRealVar(Form("nbkg_"+procname+"cat%d",catindex),Form("number of background events in signalRange "+procname+"cat%d", catindex),500,-RooNumber::infinity(),RooNumber::infinity()) ;
  RooExtendPdf* bkgpdf = new RooExtendPdf(Form("extpdf_"+procname+"cat%d_bkg",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*mdimupdf,*nbkg,"RF") ;

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
  datablind->plotOn(dtframe,Binning(tbins),Range("R1,R2"),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
  bkgpdf->plotOn(dtframe,NormRange("R1,R2"),Range("RF"));

  bkgpdf->paramOn(dtframe,Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(2))) ;
  dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextColor(kBlack);
  dtframe->getAttFill(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetFillStyle(0);
  dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextSize(0.03);
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
  if(HiggsMuMuFit::isbkg.Contains("bkg")){
      HiggsMuMuFit::bkgfit(); //bkgfit();
  }
  else if(HiggsMuMuFit::isbkg.Contains("sb")){
     HiggsMuMuFit::sbfit(); 
  }
  else HiggsMuMuFit::sigfit();
}
