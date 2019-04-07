/////////////////////////////////////////////////////////
//// 
//// Original Author : Nan Lu
////                   Caltech
//// Date Created    : Feb, 2019
/////////////////////////////////////////////////////////

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
  TString dopdf = argv[5];
  TString DNN_cut = argv[6];

  TString mass_cut = "Higgs_mass > 110 && Higgs_mass< 150 &&";
 
  TString proc = isbkgstr;

  TString sel0 = mass_cut + " disc_advNN > "+DNN_cut; //use disc_simpleNN or disc_advNN
  TString sel1 = mass_cut + " disc_advNN < "+DNN_cut+" && disc_advNN>0.0"; 
  
  vector<TString> selection = {sel0,sel1};

  int category = selection.size();
  
  HiggsMuMuFit hmm(inputFileList, outFileName, isbkgstr, funcname, selection, proc, category);
  Utils::setDefaultMinimize();
  if(dopdf.EqualTo("T")) hmm.getpdf();
  else hmm.getData(false);

  return 0;
}

void HiggsMuMuFit::sigfit()
{ 

  cout <<"sigfit: sum of up to three Gaussian functions"<<endl;
  //workspace to save pdfs, dataset, aditional variables
  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;
  //declare obs and weight vars
  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",120,110,150) ;
  mdimu->setBins(800);
  mdimu->setRange("RF",110,150);
  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-RooNumber::infinity(),RooNumber::infinity()) ;

  for(int catindex=0; catindex<ncat; catindex++){

    //load and save MC events into a RooDataSet
    RooDataSet* data = getDataSet(catindex,mdimu,evWeight,false);

    //get signal pdf
    RooAbsPdf* dimupdf = SGauss3(mdimu,catindex,data);
    if(dimupdf==NULL) dimupdf = SGauss2(mdimu,catindex,data);
    if(dimupdf==NULL){
        cout <<Form(procname+"cat%d",catindex)<<" failed, continue to next category. "<<endl;
        continue;
    }

    //expected events from MC counting events
    double sum = data->sumEntries();
    cout <<"nsig_sum: "<<sum<<endl;
    if(procname.Contains("ggH")){
      sum = sum*1.100589035;
      cout <<"scale to N3LO xs"<<endl;
    }
    cout <<"nsig_sum: "<<sum<<endl;
    RooRealVar* nsig_sum = new RooRealVar(Form(procname+"cat%d_nsig_sum",catindex),Form("number of signal events in signalRange "+procname+"cat%d", catindex),sum) ;
    nsig_sum->setConstant(true);

    //pdf x event number from MC, import this ws
    RooExtendPdf* sigpdf = new RooExtendPdf(Form(procname+"cat%d_sig_extpdf",catindex),Form("extended signal p.d.f "+procname+"cat%d", catindex),*dimupdf,*nsig_sum,"RF") ;

    //Cal FWHM
    double FWHMval = getFWHM(mdimu, dimupdf);
    cout <<"FWHM: " <<FWHMval<<endl;
    RooRealVar* FWHM = new RooRealVar(Form(procname+"cat%d_FWHM",catindex),Form("FWHM in signalRange "+procname+"cat%d", catindex),FWHMval) ;
    FWHM->setConstant(true);

    //make plot (pdf and data)
    makeplot(mdimu, dimupdf,catindex,data);
    
    //add infor into ws
    w->import(*FWHM);
    w->import(*sigpdf);
    w->import(*data);
 
    delete data;
    delete sigpdf; 
  }
  w->Print() ;
  w->writeToFile(outfile+"/pdfs/Hmm.input"+procname+"sig_13TeV.root");
}   

void HiggsMuMuFit::bias(){
  for(int i=0; i<ncat; i++) bias(i);
}

void HiggsMuMuFit::bias( int icat){
  TFile *ggH_sf = TFile::Open("output_March12_2016/pdfs/Hmm.inputggHsig_13TeV.root");
  RooWorkspace *w_ggHsig = (RooWorkspace*)ggH_sf->Get("w_all"); 
  RooExtendPdf* ggHpdf = (RooExtendPdf*)w_ggHsig->pdf(Form("ggHcat%d_sig_extpdf",icat));

  TFile *VBFH_sf = TFile::Open("output_March12_2016/pdfs/Hmm.inputVBFHsig_13TeV.root");
  RooWorkspace *w_VBFHsig = (RooWorkspace*)VBFH_sf->Get("w_all");
  RooExtendPdf* VBFHpdf = (RooExtendPdf*)w_VBFHsig->pdf(Form("VBFHcat%d_sig_extpdf",icat));

  TFile *WH_sf = TFile::Open("output_March12_2016/pdfs/Hmm.inputWHsig_13TeV.root");
  RooWorkspace *w_WHsig = (RooWorkspace*)WH_sf->Get("w_all");
  RooExtendPdf* WHpdf = (RooExtendPdf*)w_WHsig->pdf(Form("WHcat%d_sig_extpdf",icat));

  TFile *ZH_sf = TFile::Open("output_March12_2016/pdfs/Hmm.inputZHsig_13TeV.root");
  RooWorkspace *w_ZHsig = (RooWorkspace*)ZH_sf->Get("w_all");
  RooExtendPdf* ZHpdf = (RooExtendPdf*)w_ZHsig->pdf(Form("ZHcat%d_sig_extpdf",icat));
 
  cout <<"bias study for cat"<<icat<<endl;

  TFile *bkgf = TFile::Open("output_March12_2016/pdfs/Hmm.inputbkgbkg_13TeV.root");
  RooWorkspace *w_bkg = (RooWorkspace*)bkgf->Get("w_all");
  RooExtendPdf* bkgpdf = (RooExtendPdf*)w_bkg->pdf(Form("bkgcat%d_bkg_extpdf",icat));

  RooRealVar* mu = new RooRealVar(Form("mu_"+procname+"cat%d",icat),"signal strength",1.,0.0,1000.0);
  RooAbsPdf* model = new RooAddPdf(Form("sbmodel_"+procname+"cat%d", icat),"signal + bkg",RooArgList(*ggHpdf, *VBFHpdf, *WHpdf, *ZHpdf, *bkgpdf));

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
  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-100000,100000) ;

  TChain* data_chain = loader(infile.c_str(), "cattree");
  TChain* data_chain_sig = loader("/eos/cms/store/user/nlu/Hmm/VBF/categorization/VBFHToMuMu_2018_NNscore.root", "cattree");
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mdimu);
  //if(usew.EqualTo("T")) 
  obsAndWeight.add(*evWeight);

  for(int catindex=0; catindex<ncat; catindex++){
        TString cut = cuts[catindex];
        cout <<"cut: "<<cut<<endl;
        //RooCmdArg wgt_arg = usew ? RooFit::WeightVar("evt_weight") : RooCmdArg::none() ;
        RooCmdArg wgt_arg = RooFit::WeightVar("evt_weight");
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
        int status_fix = Utils::minimizeMinosTest(nll_fix,mu,1.0); 
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


void HiggsMuMuFit::bkgfit(bool isblind, bool doyield)
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

        RooDataSet* data = getDataSet(catindex,mdimu,evWeight,isblind);

        TString fitrange = "RF";
        if(isblind){
           fitrange = "R1,R2";
           cout << "DEBUG: result of fit on sideband data " << endl ;
        }
        else cout << "DEBUG: result of fit on full data " << endl ;

        RooNLLVar* nll_fix = NULL;
        if(doyield) nll_fix = (RooNLLVar*)bkgpdf->createNLL(*data, Extended(kTRUE), Range(fitrange));
        else nll_fix = (RooNLLVar*)mdimupdf->createNLL(*data);
        int status_fix = Utils::minimizeTest(nll_fix,1.0);
        if(status_fix!=0){
          cout <<"fit to data does not converge"<<endl;
          status_fix = Utils::minimizeTest(nll_fix,1.0);
          if(status_fix!=0){
             cout <<"fit to data does not converge after a second trial, try again"<<endl;
             status_fix = Utils::minimizeTest(nll_fix,1.0);
             if(status_fix!=0){
                 cout <<"fit to data does not converge after a third trial"<<endl;
             }
          } 
        }
        if(doyield) cout <<"1. expected events: "<<bkgpdf->expectedEvents(*mdimu)<<endl;

        RooBinning tbins(110,150);
        tbins.addUniform(80,110,150) ;
        RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
        data->plotOn(dtframe,Binning(tbins),Range(fitrange),Layout(0.60,0.65,0.50), Format("NEU",AutoPrecision(2)));
        if(doyield){
             bkgpdf->plotOn(dtframe,NormRange(fitrange),Range("RF"));
             bkgpdf->paramOn(dtframe,Layout(0.4, 0.88, 0.88), Format("NEU",AutoPrecision(2))) ;
             dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextColor(kBlack);
             dtframe->getAttFill(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetFillStyle(0);
             dtframe->getAttText(Form("extpdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextSize(0.03);
             cout <<"2. expected events: "<<bkgpdf->expectedEvents(*mdimu)<<endl; 
        }
        else{
             mdimupdf->plotOn(dtframe,NormRange(fitrange),Range("RF"));
             mdimupdf->paramOn(dtframe,Layout(0.4, 0.88, 0.88), Format("NEU",AutoPrecision(2))) ;
             dtframe->getAttText(Form("pdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextColor(kBlack);
             dtframe->getAttFill(Form("pdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetFillStyle(0);
             dtframe->getAttText(Form("pdf_"+procname+"cat%d_bkg_paramBox",catindex))->SetTextSize(0.03);
             nbkg->setVal(data->sumEntries()); 
             nbkg->setConstant();
             cout <<"expected events from MC: "<<nbkg->getVal()<<endl;
        }

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
      HiggsMuMuFit::bkgfit(false,false); //doblind fit? do extended fit?
  }
  else if(HiggsMuMuFit::isbkg.Contains("sb")){
     HiggsMuMuFit::sbfit(); 
  }
  else if(HiggsMuMuFit::isbkg.Contains("bias")){
     HiggsMuMuFit::bias();
  }
  else HiggsMuMuFit::sigfit();
}

void HiggsMuMuFit::getData(bool isblind){
  //workspace to save dataset
  RooWorkspace *w = new RooWorkspace("w_all","workspace") ;
  //declare obs and weight vars
  RooRealVar* mdimu = new RooRealVar("Higgs_mass","Higgs_mass",120,110,150) ;
  mdimu->setBins(4000);
  mdimu->setRange("RF",110,150);
  RooRealVar* evWeight = new RooRealVar("evt_weight","evt_weight",1,-RooNumber::infinity(),RooNumber::infinity()) ;

  for(int catindex=0; catindex<ncat; catindex++){
     //load and save MC events into a RooDataSet
     RooDataSet* data = getDataSet(catindex,mdimu,evWeight,isblind);
     RooDataHist* databinned = data->binnedClone();
     cout <<catindex<<" data number from unbinned dataset: "<<data->sumEntries()<<endl;
     cout <<catindex<<" data infor from unbinned dataset"<<endl;
     databinned->Print("v");
     w->import(*data);
     w->import(*databinned);
  }
  w->Print() ;
  w->writeToFile(outfile+"/dataset/Hmm.input"+procname+"_13TeV.root");
}
