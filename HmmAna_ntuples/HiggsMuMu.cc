#define HiggsMuMu_cxx
#include "HiggsMuMu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>
int main(int argc, char* argv[])
{

  if(argc < 3) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  HiggsMuMu hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << endl;
  hmm.EventLoop(data, isData);

  return 0;
}

void HiggsMuMu::EventLoop(const char *data,const char *isData)
{  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double evt_wt;
      if(*isData=='T'){evt_wt=1.;}
      else evt_wt=t_genWeight;
      if(t_diMuon_mass<120. || t_diMuon_mass>130.)h_diMuon_mass->Fill(t_diMuon_mass,evt_wt);
      if((t_diMuon_mass>110. && t_diMuon_mass<120.) || (t_diMuon_mass>130. && t_diMuon_mass<150.)){
	h_mu1pt->Fill((*t_Mu_pt)[0],evt_wt);
	h_mu2pt->Fill((*t_Mu_pt)[1],evt_wt);
	h_mu1eta->Fill((*t_Mu_eta)[0],evt_wt);
	h_mu2eta->Fill((*t_Mu_eta)[1],evt_wt);
	h_mu1phi->Fill((*t_Mu_phi)[0],evt_wt);
	h_mu2phi->Fill((*t_Mu_phi)[1],evt_wt);
	double dR= DeltaR((*t_Mu_eta)[0],(*t_Mu_phi)[0],(*t_Mu_eta)[1],(*t_Mu_phi)[1]);
	h_mu1mu2dR->Fill(dR,evt_wt);
	double dPhi= DeltaPhi((*t_Mu_phi)[1],(*t_Mu_phi)[0]);
	h_mu1mu2dPhi->Fill(dPhi,evt_wt);
	// if (Cut(ientry) < 0) continue;
	
	h_diMuon_pt->Fill(t_diMuon_pt,evt_wt);
	h_diMuon_eta->Fill(t_diMuon_eta,evt_wt);
	h_diMuon_phi->Fill(t_diMuon_phi,evt_wt);
	h_diMuon_mass_110To150->Fill(t_diMuon_mass,evt_wt);
	if(t_diMuon_mass>110.&& t_diMuon_mass<120.)h_diMuon_mass_110To120->Fill(t_diMuon_mass,evt_wt);
	if(t_diMuon_mass>130. && t_diMuon_mass<150.)h_diMuon_mass_130To150->Fill(t_diMuon_mass,evt_wt);
	
	h_MET_pt->Fill(t_MET_pt);
	h_METphi->Fill(t_MET_phi);
	h_MET_sumEt->Fill(t_MET_sumEt);

	h_Njet->Fill(t_nJet);
	if(t_nJet>0){
	  h_j1pt->Fill((*t_Jet_pt)[0],evt_wt);
	  h_j1eta->Fill((*t_Jet_eta)[0],evt_wt);
	  h_j1phi->Fill((*t_Jet_phi)[0],evt_wt);
	}

	if(t_nJet>1){
	  h_j2pt->Fill((*t_Jet_pt)[1],evt_wt);
	  h_j2eta->Fill((*t_Jet_eta)[1],evt_wt);
	  h_j2phi->Fill((*t_Jet_phi)[1],evt_wt);

	  dR= DeltaR((*t_Jet_eta)[0],(*t_Jet_phi)[0],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
	  h_j1j2dR->Fill(dR,evt_wt);
	  dPhi= DeltaPhi((*t_Jet_phi)[0],(*t_Jet_phi)[1]);
	  h_j1j2dPhi->Fill(dPhi,evt_wt);

	  h_dijet_pt->Fill(t_diJet_pt,evt_wt);
	  h_dijet_eta->Fill(t_diJet_eta,evt_wt);
	  h_dijet_phi->Fill(t_diJet_phi,evt_wt);
	  h_Mjj->Fill(t_diJet_mass,evt_wt);
	}
      }
   }
}

