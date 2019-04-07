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
#include "Utils.h"

using namespace RooFit ;
using namespace std ;

class HiggsMuMuFit{
public :

  HiggsMuMuFit(TString inputFileList="foo.txt", TString outFileName="histo.root", TString isbkgstr="F", TString funcname="test", vector<TString> selection={"cat_index==7"}, TString proc="ggH", int category=1);
  virtual ~HiggsMuMuFit();
  void     sigfit();
  void     sigfit_DCB();
  void     bkgfit(bool isblind, bool doyield);
  void     sbfit();
  RooAbsPdf* BWZRedux(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZReduxB4(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZReduxB3(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZGamma(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZ(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZPol1(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* BWZPol2(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* SExp(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex);
  RooAbsPdf* SGauss3(RooRealVar* mdimu,int catindex, RooDataSet* dataset);
  RooAbsPdf* SGauss2(RooRealVar* mdimu,int catindex, RooDataSet* dataset);
  RooAbsPdf* SGauss(RooRealVar* mdimu,int catindex, double m1, double m2, double s1, double s2, double f1);
  RooAbsPdf* bkgfunc(RooRealVar*mdimu, RooFormulaVar* mdimu_range, int catindex, TString func_name);
  double getFWHM(RooRealVar *mass, RooAbsPdf *pdf, double wmin, double wmax, double step);
  void makeplot(RooRealVar* mdimu, RooAbsPdf* dimupdf, int catindex, RooDataSet* dataset);
  RooDataSet* getDataSet(int catindex, RooRealVar* mdimu, RooRealVar* evWeight, bool isblind);
  void     getData(bool isblind);
  void     bias(); 
  void     bias(int icat);
  void     getpdf();
  void     set(TString inputFileList, TString outFileName, TString isbkgstr, TString funcname, vector<TString> selection, TString proc, int category);
  TChain* loader(const string& inFile_name, const string& chain_name);
  void InitTreeVars();
  void BookTreeBranches();

  //TFile* oFile;
  TTree* outtree;
  int t_category;
  double t_mu;
  double t_errlo;
  double t_errhi;
  double t_err;
  double t_nbkg;
  double t_nsig;

private:
  string infile;
  TString outfile;
  vector<TString> cuts;
  TString isbkg;
  vector<TString> funclist;
  TString procname;
  int ncat;
};

HiggsMuMuFit::HiggsMuMuFit(TString inputFileList, TString outFileName, TString isbkgstr, TString funcname, vector<TString> selection, TString proc, int category)
{
  infile = inputFileList.Data();
  outfile = outFileName;
  cuts = selection;
  isbkg = isbkgstr;
  funclist=Utils::SplitString(funcname,',');
  procname = proc;
  ncat = category; 
  //oFile = new TFile(outFileName+"/"+proc+"_result.root", "recreate");
  InitTreeVars();
  BookTreeBranches();
}

void HiggsMuMuFit::set(TString inputFileList, TString outFileName, TString isbkgstr, TString funcname, vector<TString> selection, TString proc, int category)
{
  infile = inputFileList;
  outfile = outFileName;
  cuts = selection;
  isbkg = isbkgstr;
  funclist=Utils::SplitString(funcname,',');
  procname = proc;
  ncat = category;
}

void HiggsMuMuFit::InitTreeVars(){
  t_category = -999;
  t_mu = -999.;
  t_errhi = -999.;
  t_errlo = -999.;
  t_err = -999.;
  t_nbkg = -999.;
  t_nsig = -999.;
}

void HiggsMuMuFit::BookTreeBranches(){
  outtree = new TTree("result","result");
  outtree->Branch("category", &t_category,"category/i");
  outtree->Branch("mu", &t_mu,"mu/D");
  outtree->Branch("errhi", &t_errhi,"errhi/D");
  outtree->Branch("errlo", &t_errlo,"errlo/D");
  outtree->Branch("err", &t_err,"err/D");
  outtree->Branch("nbkg", &t_nbkg,"nbkg/D");
  outtree->Branch("nsig", &t_nsig,"nsig/D");
}

TChain* HiggsMuMuFit::loader(const string& inFile_name, const string& chain_name)
{
    TChain* chain = new TChain(chain_name.c_str());
    TString in_name(inFile_name);
    if(in_name.Contains("root")) {
        chain->Add(inFile_name.c_str());
        cout << "single root file, total events: " << chain->GetEntries() << " in " << inFile_name.c_str() << endl;
        return chain;
    }
    fstream input(inFile_name.c_str(), fstream::in);
    string file_name;
    int ncounter = 0;
    while (input >> file_name){
        cout << "adding: " << file_name << endl;
        chain->Add(file_name.c_str());
        ncounter ++;
    }
    cout << "total events: " << chain->GetEntries() << " in " << ncounter << " files." << endl;
    input.close();
    return chain;
}

HiggsMuMuFit::~HiggsMuMuFit()
{

/*
oFile->cd();
outtree->Write("", TObject::kOverwrite);
oFile->Write();
oFile->Close();
*/
}

RooAbsPdf* HiggsMuMuFit::BWZRedux(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-20,20) ;
    RooRealVar* btwzr_a2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),0.1,-100,100) ;
    RooRealVar* btwzr_a3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),0.1,-20,20) ;
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3));

  return mdimupdf;
}

RooAbsPdf* HiggsMuMuFit::BWZReduxB4(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
    RooRealVar* btwzr_a2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),0.1,-20,20) ;
    RooRealVar* btwzr_a3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),0.1,-20,20) ;
    RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b4 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b4_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b4_cat%d",catindex),0.01,-1.0,1.0) ;
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*5*pow(@4,4)+20*@6*pow(@4,3)*(1-@4)+@7*30*@4*@4*(1-@4)*(1-@4) + @8 *20*@4*pow(1-@4,3)+(-@5-@6-@7-@8)*5*pow(1-@4,4))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3,*b4));

  return mdimupdf;
}

RooAbsPdf* HiggsMuMuFit::BWZReduxB3(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
    RooRealVar* btwzr_a2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a2_cat%d",catindex),0.1,-20,20) ;
    RooRealVar* btwzr_a3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a3_cat%d",catindex),0.1,-20,20) ;
    RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),0.01,-1.0,1.0) ;
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@2*@0/100. +(@0/100.)*(@0/100.)*@3 )/(TMath::Power((@0-91.2),@1)+TMath::Power(2.5/2.,@1)) * (1+ @5*5*pow(@4,3)+20*@6*pow(@4,2)*(1-@4)+@7*30*@4*(1-@4)*(1-@4) + (-@5-@6-@7)*20*pow(1-@4,3))",RooArgList(*mdimu,*btwzr_a1,*btwzr_a2,*btwzr_a3,*mdimu_range,*b1,*b2,*b3));

  return mdimupdf;
}


RooAbsPdf* HiggsMuMuFit::BWZGamma(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
    RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;

    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@1*@0/100. )*2.5*@2/(TMath::Power((@0-91.2),2)+TMath::Power(2.5/2.,2)) + (1-@2)*TMath::Exp(@1*@0/100.)/(@0*@0)",RooArgList(*mdimu,*btwzr_a1,*b1));

    return mdimupdf;
}

RooAbsPdf* HiggsMuMuFit::SExp(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),5.0,-100.0,100.0);
    RooRealVar* b2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),5.0,-100.0,100.0);
    RooRealVar* b3 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b3_cat%d",catindex),5.0,-100.0,100.0);
    RooRealVar* a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_a1_cat%d",catindex),0.5,0.0,1.0);
    RooRealVar* a2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_a2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_a2_cat%d",catindex),0.5,0.0,1.0);
    //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(-@2*@0/100.)*@1+TMath::Exp(-@3*@0/100.)*(1-@1)",RooArgList(*mdimu,*b1,*b2,*b3));
    //RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"(@3)*TMath::Exp(-@1*@0/100.)+(1.0-@3)*TMath::Exp(-@2*@0/100.)",RooArgList(*mdimu,*b1,*b2,*a1));
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"(@4)*TMath::Exp(-@1*@0/100.)+(@5)*TMath::Exp(-@2*@0/100.)+(1.0-@4-@5)*TMath::Exp(-@3*@0/100.)",RooArgList(*mdimu,*b1,*b2,*b3,*a1,*a2));

  return mdimupdf;
}

RooAbsPdf* HiggsMuMuFit::BWZ(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@1*@0/100. )*2.5/(TMath::Power((@0-91.2),2)+TMath::Power(2.5/2.,2))",RooArgList(*mdimu,*btwzr_a1));

  return mdimupdf;
}


RooAbsPdf* HiggsMuMuFit::BWZPol1(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

   RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
   RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;
   RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@1*@0/100. )*2.5*(@2*@3+1-@3)/(TMath::Power((@0-91.2),2)+TMath::Power(2.5/2.,2))",RooArgList(*mdimu,*btwzr_a1,*mdimu_range,*b1));

  return mdimupdf;
}


RooAbsPdf* HiggsMuMuFit::BWZPol2(RooRealVar* mdimu, RooFormulaVar* mdimu_range, int catindex){

    RooRealVar* btwzr_a1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_btwzr_a1_cat%d",catindex),0.1,-5,5) ;
    RooRealVar* b1 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b1_cat%d",catindex),0.01,-1.0,1.0) ;
    RooRealVar* b2 = new RooRealVar(TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),TString::Format("CMS_Hmm_"+procname+"_b2_cat%d",catindex),0.01,-1.0,1.0) ;
    RooAbsPdf* mdimupdf = new RooGenericPdf(TString::Format("pdf_"+procname+"cat%d_bkg",catindex),"TMath::Exp(@1*@0/100. )*2.5*(@2*@2*@3+@2*@4+1-@3-@4)/(TMath::Power((@0-91.2),2)+TMath::Power(2.5/2.,2))",RooArgList(*mdimu,*btwzr_a1,*mdimu_range,*b1,*b2));
  return mdimupdf;
}

RooAbsPdf* HiggsMuMuFit::bkgfunc(RooRealVar*mdimu, RooFormulaVar* mdimu_range, int catindex, TString func_name){
        RooAbsPdf* bkgpdf = NULL;
        if(func_name.EqualTo("BWZRedux")){
            bkgpdf = HiggsMuMuFit::BWZRedux( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZReduxB4")){
            bkgpdf = HiggsMuMuFit::BWZReduxB4( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZReduxB3")){
            bkgpdf = HiggsMuMuFit::BWZReduxB3( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZ")){
            bkgpdf = HiggsMuMuFit::BWZ( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZPol1")){
            bkgpdf = HiggsMuMuFit::BWZPol1( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZPol2")){
            bkgpdf = HiggsMuMuFit::BWZPol2( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("BWZGamma")){
            bkgpdf = HiggsMuMuFit::BWZGamma( mdimu, mdimu_range, catindex);
        }
        else if(func_name.EqualTo("SExp")){
            bkgpdf = HiggsMuMuFit::SExp( mdimu, mdimu_range, catindex);
        }
        else{
          cout <<"function does not find, error"<<endl;
        }
        return bkgpdf;
}

RooAbsPdf* HiggsMuMuFit::SGauss(RooRealVar* mdimu,int catindex, double m1, double m2, double s1, double s2, double f1){
    RooRealVar* mean1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean1_cat%d",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex), m1, m1, m1) ;
    RooRealVar* mean2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean2_cat%d",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex), m2, m2, m2) ;
    RooRealVar* sigma1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma1_cat%d",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex), s1, s1, s1) ;
    RooRealVar* sigma2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma2_cat%d",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex), s2, s2, s2) ;
    RooGaussian* sig1 = new RooGaussian(Form("CMS_Hmm_"+procname+"_sig1_cat%d",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
    RooGaussian* sig2 = new RooGaussian(Form("CMS_Hmm_"+procname+"_sig2_cat%d",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;
    RooRealVar* sig1frac = new RooRealVar(Form("CMS_Hmm_"+procname+"_sig1frac_cat%d",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),f1,f1,f1) ;   

    mean1->setConstant(true);
    mean2->setConstant(true);
    sigma1->setConstant(true);
    sigma2->setConstant(true);
    sig1frac->setConstant(true);

    RooAddPdf* dimupdf = new RooAddPdf(Form(procname+"cat%d_sig_pdf",catindex),Form(procname+"cat%d_sig_pdf",catindex),RooArgList(*sig1,*sig2),RooArgList(*sig1frac)) ;
    return dimupdf;

}

RooAbsPdf* HiggsMuMuFit::SGauss3(RooRealVar* mdimu,int catindex, RooDataSet* data){

    RooRealVar* mean1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean1_cat%d",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
    RooRealVar* mean2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean2_cat%d",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
    RooRealVar* mean3 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean3_cat%d",catindex),Form("mean3 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;

    RooRealVar* sigma1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma1_cat%d",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex),3, 0.0, 10) ;
    RooRealVar* sigma2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma2_cat%d",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex),5, 0.0, 40) ;
    RooRealVar* sigma3 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma3_cat%d",catindex),Form("sigma3 of gaussians for "+procname+"cat%d", catindex),10, 0.0, 40) ;

    RooGaussian* sig1 = new RooGaussian(Form(procname+"cat%d_sig1",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
    RooGaussian* sig2 = new RooGaussian(Form(procname+"cat%d_sig2",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;
    RooGaussian* sig3 = new RooGaussian(Form(procname+"cat%d_sig3",catindex),Form("Signal component 3 for "+procname+"cat%d", catindex),*mdimu, *mean3, *sigma3) ;
 
    RooRealVar* sig1frac = new RooRealVar(Form(procname+"cat%d_sig1frac",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),0.8,0.,1.) ;
    RooRealVar* sig2frac = new RooRealVar(Form(procname+"cat%d_sig2frac",catindex),Form("signal fraction of component 2 for "+procname+"cat%d", catindex),0.1,0.,1.) ;
    RooAddPdf* dimupdf = new RooAddPdf(Form("pdf_"+procname+"cat%d_sig",catindex),Form("pdf_"+procname+"cat%d_sig",catindex),RooArgList(*sig1,*sig2, *sig3),RooArgList(*sig1frac,*sig2frac)) ;

    RooNLLVar* nll = (RooNLLVar*)dimupdf->createNLL(*data);
    nll->enableOffsetting(true);
    int status_fit = Utils::minimizeTest(nll,1.0);
    if((status_fit!=0) && (status_fit!=1)){
       cout <<"fit to data does not converge for 3 Gauss, try again"<<endl;
       status_fit = Utils::minimizeTest(nll,1.0);
       if((status_fit!=0) && (status_fit!=1)){
          cout <<"fit to data does not converge"<<endl;
          return NULL;
       }
    }

    sig1frac->setConstant(true);
    sig2frac->setConstant(true);
    sigma1->setConstant(true);
    sigma2->setConstant(true);
    sigma3->setConstant(true);
    mean1->setConstant(true);
    mean2->setConstant(true);
    mean3->setConstant(true);
    return dimupdf; 
}

RooAbsPdf* HiggsMuMuFit::SGauss2(RooRealVar* mdimu,int catindex, RooDataSet* data){

    RooRealVar* mean1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean1_cat%d",catindex),Form("mean1 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;
    RooRealVar* mean2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_mean2_cat%d",catindex),Form("mean2 of gaussians for "+procname+"cat%d", catindex),125, 110, 145) ;

    RooRealVar* sigma1 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma1_cat%d",catindex),Form("sigma1 of gaussians for "+procname+"cat%d", catindex),3, 0.0, 10) ;
    RooRealVar* sigma2 = new RooRealVar(Form("CMS_Hmm_"+procname+"_sigma2_cat%d",catindex),Form("sigma2 of gaussians for "+procname+"cat%d", catindex),5, 0.0, 40) ;

    RooGaussian* sig1 = new RooGaussian(Form(procname+"cat%d_sig1",catindex),Form("Signal component 1 for "+procname+"cat%d", catindex),*mdimu, *mean1, *sigma1) ;
    RooGaussian* sig2 = new RooGaussian(Form(procname+"cat%d_sig2",catindex),Form("Signal component 2 for "+procname+"cat%d", catindex),*mdimu, *mean2, *sigma2) ;

    RooRealVar* sig1frac = new RooRealVar(Form(procname+"cat%d_sig1frac",catindex),Form("signal fraction of component 1 for "+procname+"cat%d", catindex),0.8,0.,1.) ;
    RooAddPdf* dimupdf = new RooAddPdf(Form("pdf_"+procname+"cat%d_sig",catindex),Form("pdf_"+procname+"cat%d_sig",catindex),RooArgList(*sig1,*sig2),RooArgList(*sig1frac)) ;

    RooNLLVar* nll = (RooNLLVar*)dimupdf->createNLL(*data);
    nll->enableOffsetting(true);
    int status_fit = Utils::minimizeTest(nll,1.0);
    if(status_fit!=0){
       cout <<"fit to data does not converge for 2 Gauss, try again"<<endl;
       status_fit = Utils::minimizeTest(nll,1.0);
       if(status_fit!=0){
          cout <<"fit to data does not converge"<<endl;
          return NULL;
       }
    }

    sig1frac->setConstant(true);
    sigma1->setConstant(true);
    sigma2->setConstant(true);
    mean1->setConstant(true);
    mean2->setConstant(true);

    return dimupdf;
}

double HiggsMuMuFit::getFWHM(RooRealVar *mass, RooAbsPdf *pdf, double wmin=110., double wmax=130., double step=0.025) {
//vector<double> getFWHM(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, double wmin=110., double wmax=130., double step=0.025) {

  cout << "Computing FWHM...." << endl;
  double nbins = (wmax-wmin)/step;
  TH1F *h = new TH1F("h","h",int(floor(nbins+0.5)),wmin,wmax);
  //if (data){
  //  pdf->fillHistogram(h,RooArgList(*mass),data->sumEntries());
  //}
  //else {
    pdf->fillHistogram(h,RooArgList(*mass));
  //}

  double hm = h->GetMaximum()*0.5;
  double low = h->GetBinCenter(h->FindFirstBinAbove(hm));
  double high = h->GetBinCenter(h->FindLastBinAbove(hm));

  double val = high-low;
  cout << "FWHM: [" << low << "-" << high << " = " << high-low <<"] Max = " << hm << endl;
  /*
  vector<double> result;
  result.push_back(low);
  result.push_back(high);
  result.push_back(hm);
  result.push_back(h->GetBinWidth(1));
  */

  delete h;
  //return result;
  return val;
}

void HiggsMuMuFit::makeplot(RooRealVar* mdimu, RooAbsPdf* dimupdf,int catindex, RooDataSet* data){

    RooBinning tbins(110,150);
    tbins.addUniform(80,110,150) ;
    RooPlot* dtframe = mdimu->frame(Range(110,150),Title("m(#mu#mu) distribution"));
    cout <<"p1"<<endl;
    data->plotOn(dtframe,Binning(tbins));
    cout <<"p2"<<endl;
    RooArgSet* check = dimupdf->getComponents();
    check->Print();
    TIterator* iter = check->createIterator();
    TObject* arg = NULL;
    vector<string> nComp;
    nComp.clear();
    int icolor = 1;
    while((arg=(TObject*)iter->Next())){
        string str(arg->GetName());
        icolor++;
        cout <<"plot: "<<str<<endl;
        nComp.push_back(str);
        if(icolor==2) dimupdf->plotOn(dtframe,Components(str.data()),LineColor(icolor),Name(str.data()));
        else{
          if(icolor==10) dimupdf->plotOn(dtframe,Components(str.data()),LineColor(kOrange+7),Name(str.data()),LineStyle(kDashed));
          else dimupdf->plotOn(dtframe,Components(str.data()),LineColor(icolor),Name(str.data()),LineStyle(kDashed));
        }
    }
    /*
 *     args = RooArgSet (parameters...);
 *         dimupdf->paramOn(dtframe,RooFit.Parameters(args), Layout(0.55, 0.88, 0.88), Format("NEU",AutoPrecision(2)));
 *             dtframe->getAttText(Form("pdf_"+procname+"cat%d_sig_paramBox",catindex))->SetTextColor(kBlack);
 *                 dtframe->getAttFill(Form("pdf_"+procname+"cat%d_sig_paramBox",catindex))->SetFillStyle(0);
 *                     dtframe->getAttText(Form("pdf_"+procname+"cat%d_sig_paramBox",catindex))->SetTextSize(0.03);
 *                         */

    TCanvas* c1 = new TCanvas("c1","c1");
    dtframe->Draw();
    c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m.png",catindex)));
    gPad->SetLogy();
    double maxv = dtframe->GetMaximum();
    dtframe->SetMinimum(0.01*maxv);
    dtframe->Draw();
    c1->SaveAs(TString::Format(Form(outfile+"/plots/Hmm."+procname+"sig_13TeV.cat%d_m_log.png",catindex)));

}

RooDataSet* HiggsMuMuFit::getDataSet(int catindex, RooRealVar* mdimu, RooRealVar* evWeight, bool isblind){
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mdimu);
    obsAndWeight.add(*evWeight);

    TChain* data_chain = loader(infile.c_str(), "cattree");
    TString cut = cuts[catindex];
    if(isblind) cut = cut + "&& (Higgs_mass>130 || Higgs_mass<120)";
    cout <<"cut: "<<cut<<endl;
    TTree* cutChain = data_chain->CopyTree(cut);
    string dataset_ch_name(Form("Sig_"+procname+"_cat%d",catindex));
    RooCmdArg wgt_arg = RooFit::WeightVar("evt_weight");
    RooDataSet* data = new RooDataSet(dataset_ch_name.c_str(), dataset_ch_name.c_str(), RooArgSet(obsAndWeight), RooFit::Import(*cutChain), wgt_arg);
    delete data_chain;
    delete cutChain;
    return data;
}

