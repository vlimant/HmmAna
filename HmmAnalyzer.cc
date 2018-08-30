#define HmmAnalyzer_cxx
#include "HmmAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

int main(int argc, char* argv[])
{

  if(argc < 4) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  HmmAnalyzer Hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << endl;
  Hmm.EventLoop(data, isData);

  return 0;
}

void HmmAnalyzer::EventLoop(const char *data,const char *isData)
{ 
  if (fChain == 0) return;
  //clearTreeVectors();
  //cout<<"cleared tree vectors\n";
  //BookTreeBranches();
  //cout<<"booked tree branches\n";
  Long64_t nentries = fChain->GetEntriesFast();
  
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      bool trig_decision = false;

      if( HLT_IsoMu27==1 /* || HLT_IsoTkMu27_v*==1*/) trig_decision =true;

      bool goodLumi= false;
      if (isData=="T"){//figure out good lumi from json
      }
      else  goodLumi = true;
      
      bool run_muChecks =false; 
      if(Flag_METFilters && trig_decision && goodLumi && PV_ndof>4 && fabs(PV_z)<24. && PV_npvsGood>0) run_muChecks =true;

      bool Event_sel= false;
      if(run_muChecks){
	if(Muon_charge[0]*Muon_charge[1]== -1 && Muon_pt[0]>30.) Event_sel =true;
      }
      if(Event_sel){

	t_run =run;
	t_luminosityBlock=luminosityBlock;
	t_event=event;
	

	for(int i=0;i<nMuon;i++){
	  if(Muon_pt[i]>20. && fabs(Muon_eta[i])<2.4 && Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.25){
	    t_Mu_charge->push_back(Muon_charge[i]);   
	    t_Mu_pt->push_back(Muon_pt[i]);   
	    t_Mu_ptErr->push_back(Muon_ptErr[i]);   
	    t_Mu_phi->push_back(Muon_phi[i]);   
	    t_Mu_eta->push_back(Muon_eta[i]);   
	    t_Mu_mass->push_back(Muon_mass[i]);  
	    t_Mu_dxy->push_back(Muon_dxy[i]);   
	    t_Mu_dxyErr->push_back(Muon_dxyErr[i]);   
	    t_Mu_dz->push_back(Muon_dz[i]);   
	    t_Mu_dzErr->push_back(Muon_dzErr[i]); 
	    t_Mu_pfRelIso03_all->push_back(Muon_pfRelIso03_all[i]);   
	    t_Mu_pfRelIso03_chg->push_back(Muon_pfRelIso03_chg[i]);   
	    t_Mu_pfRelIso04_all->push_back(Muon_pfRelIso04_all[i]);   
	    t_Mu_tightCharge->push_back(Muon_tightCharge[i]);   
	    t_Mu_isPFcand->push_back(Muon_isPFcand[i]);   
	    t_Mu_mediumId->push_back(Muon_mediumId[i]);   
	    t_Mu_softId->push_back(Muon_softId[i]);   
	    t_Mu_tightId->push_back(Muon_tightId[i]);    
	    t_Mu_nStations->push_back(Muon_nStations[i]);   
	    t_Mu_nTrackerLayers->push_back(Muon_nTrackerLayers[i]);   
	  }
	}

	if(Muon_pt[1]>20. && fabs(Muon_eta[1])<2.4 && Muon_mediumId[1] && Muon_pfRelIso04_all[1] < 0.25){
	  if(Muon_pt[0]>30. && fabs(Muon_eta[0])<2.4 && Muon_mediumId[0] && Muon_pfRelIso04_all[0] < 0.25){
	    TLorentzVector dimu, mu1,mu2;
	    mu1.SetPtEtaPhiM(Muon_pt[0],Muon_eta[0],Muon_phi[0],Muon_mass[0]);
	    mu2.SetPtEtaPhiM(Muon_pt[1],Muon_eta[1],Muon_phi[1],Muon_mass[1]);
	    dimu=mu1+mu2;
	    t_diMuon_pt = dimu.Pt();
	    t_diMuon_eta= dimu.Eta();
	    t_diMuon_phi= dimu.Phi();
	    t_diMuon_mass= dimu.M();
	  }
	}
	for (int j =0;j<nJet;j++){
	  if(Jet_muonIdx1[j]==-1 && Jet_muonIdx2[j]==-1 && Jet_pt[j]>30. && fabs(Jet_eta[j])<4.7 && Jet_jetId[j]>=2 && Jet_puId[j]>=1){
	    t_Jet_area->push_back(Jet_area[j]);
	    t_Jet_btagCMVA->push_back(Jet_btagCMVA[j]);   
	    t_Jet_btagCSVV2->push_back(Jet_btagCSVV2[j]);   
	    t_Jet_btagDeepB->push_back(Jet_btagDeepB[j]);   
	    t_Jet_btagDeepC->push_back(Jet_btagDeepC[j]);   
	    t_Jet_btagDeepFlavB->push_back(Jet_btagDeepFlavB[j]);   
	    t_Jet_chEmEF->push_back(Jet_chEmEF[j]);   
	    t_Jet_chHEF->push_back(Jet_chHEF[j]);   
	    t_Jet_eta->push_back(Jet_eta[j]);   
	    t_Jet_mass->push_back(Jet_mass[j]);   
	    t_Jet_neEmEF->push_back(Jet_neEmEF[j]);   
	    t_Jet_neHEF->push_back(Jet_neHEF[j]);   
	    t_Jet_phi->push_back(Jet_phi[j]);   
	    t_Jet_pt->push_back(Jet_pt[j]);   
	    t_Jet_qgl->push_back(Jet_qgl[j]);   
	    t_Jet_jetId->push_back(Jet_jetId[j]);   
	    t_Jet_nConstituents->push_back(Jet_nConstituents[j]);   
	    t_Jet_nElectrons->push_back(Jet_nElectrons[j]);   
	    t_Jet_nMuons->push_back(Jet_nMuons[j]);   
	    t_Jet_puId->push_back(Jet_puId[j]);   

	    if(Jet_btagDeepB[j]>0.4941){
	      t_bJet_area->push_back(Jet_area[j]);
	      t_bJet_btagCMVA->push_back(Jet_btagCMVA[j]);   
	      t_bJet_btagCSVV2->push_back(Jet_btagCSVV2[j]);   
	      t_bJet_btagDeepB->push_back(Jet_btagDeepB[j]);   
	      t_bJet_btagDeepC->push_back(Jet_btagDeepC[j]);   
	      t_bJet_btagDeepFlavB->push_back(Jet_btagDeepFlavB[j]);   
	      t_bJet_chEmEF->push_back(Jet_chEmEF[j]);   
	      t_bJet_chHEF->push_back(Jet_chHEF[j]);   
	      t_bJet_eta->push_back(Jet_eta[j]);   
	      t_bJet_mass->push_back(Jet_mass[j]);   
	      t_bJet_neEmEF->push_back(Jet_neEmEF[j]);   
	      t_bJet_neHEF->push_back(Jet_neHEF[j]);   
	      t_bJet_phi->push_back(Jet_phi[j]);   
	      t_bJet_pt->push_back(Jet_pt[j]);   
	      t_bJet_qgl->push_back(Jet_qgl[j]);   
	      t_bJet_jetId->push_back(Jet_jetId[j]);   
	      t_bJet_nConstituents->push_back(Jet_nConstituents[j]);   
	      t_bJet_nElectrons->push_back(Jet_nElectrons[j]);   
	      t_bJet_nMuons->push_back(Jet_nMuons[j]);   
	      t_bJet_puId->push_back(Jet_puId[j]);   

	  
	    }
 
	  }
	}

	for(int i=0;i<nElectron;i++){
	  t_El_charge->push_back(Electron_charge[i]);
	  t_El_pt->push_back(Electron_pt[i]);
	  t_El_phi->push_back(Electron_phi[i]);
	  t_El_eta->push_back(Electron_eta[i]);   
	  t_El_mass->push_back(Electron_mass[i]);
	  t_El_cutBased->push_back(Electron_cutBased[i]);   
	  t_El_tightCharge->push_back(Electron_tightCharge[i]);   
	  t_El_cutBased_HEEP->push_back(Electron_cutBased_HEEP[i]);   
	  t_El_isPFcand->push_back(Electron_isPFcand[i]);   
	  t_El_pfRelIso03_all->push_back(Electron_pfRelIso03_all[i]);   
	  t_El_pfRelIso03_chg->push_back(Electron_pfRelIso03_chg[i]);      
	  t_El_dxy->push_back(Electron_dxy[i]);   
	  t_El_dxyErr->push_back(Electron_dxyErr[i]);   
	  t_El_dz->push_back(Electron_dz[i]);   
	  t_El_dzErr->push_back(Electron_dzErr[i]);   
	  t_Electron_mvaFall17Iso_WP80->push_back(Electron_mvaFall17Iso_WP80[i]);
	  t_Electron_mvaFall17Iso_WP90->push_back(Electron_mvaFall17Iso_WP90[i]);
	  t_Electron_mvaFall17Iso_WPL->push_back(Electron_mvaFall17Iso_WPL[i]);
	  t_Electron_mvaFall17noIso_WP80->push_back(Electron_mvaFall17noIso_WP80[i]);
	  t_Electron_mvaFall17noIso_WP90->push_back(Electron_mvaFall17noIso_WP90[i]);
	  t_Electron_mvaFall17noIso_WPL->push_back(Electron_mvaFall17noIso_WPL[i]);
  
	}

	t_MET_pt = MET_pt;
	t_MET_phi = MET_phi;
	t_MET_sumEt  = MET_sumEt;


	for(int i=0;i<nFatJet;i++){
	  t_FatJet_area->push_back(FatJet_area[i]);  
	  t_FatJet_btagCMVA->push_back(FatJet_btagCMVA[i]);  
	  t_FatJet_btagCSVV2->push_back(FatJet_btagCSVV2[i]);  
	  t_FatJet_btagDeepB->push_back(FatJet_btagDeepB[i]);  
	  t_FatJet_eta->push_back(FatJet_eta[i]);  
	  t_FatJet_mass->push_back(FatJet_mass[i]);  
	  t_FatJet_msoftdrop->push_back(FatJet_msoftdrop[i]);  
	  t_FatJet_n2b1->push_back(FatJet_n2b1[i]);  
	  t_FatJet_n3b1->push_back(FatJet_n3b1[i]);  
	  t_FatJet_phi->push_back(FatJet_phi[i]);  
	  t_FatJet_pt->push_back(FatJet_pt[i]);  
	  t_FatJet_tau1->push_back(FatJet_tau1[i]);  
	  t_FatJet_tau2->push_back(FatJet_tau2[i]);  
	  t_FatJet_tau3->push_back(FatJet_tau3[i]);  
	  t_FatJet_tau4->push_back(FatJet_tau4[i]);  
	  t_FatJet_jetId->push_back(FatJet_tau4[i]);  
	  t_FatJet_subJetIdx1->push_back(FatJet_subJetIdx1[i]);  
	  t_FatJet_subJetIdx2->push_back(FatJet_subJetIdx2[i]);  
	}
	for(int i=0;i<nSubJet;i++){
	  t_SubJet_btagCMVA->push_back(SubJet_btagCMVA[i]);   
	  t_SubJet_btagCSVV2->push_back(SubJet_btagCSVV2[i]);   
	  t_SubJet_btagDeepB->push_back(SubJet_btagDeepB[i]);   
	  t_SubJet_eta->push_back(SubJet_eta[i]);   
	  t_SubJet_mass->push_back(SubJet_mass[i]);   
	  t_SubJet_n2b1->push_back(SubJet_n2b1[i]);   
	  t_SubJet_n3b1->push_back(SubJet_n3b1[i]);   
	  t_SubJet_phi->push_back(SubJet_phi[i]);   
	  t_SubJet_pt->push_back(SubJet_pt[i]);   
	  t_SubJet_tau1->push_back(SubJet_tau1[i]);   
	  t_SubJet_tau2->push_back(SubJet_tau2[i]);   
	  t_SubJet_tau3->push_back(SubJet_tau3[i]);   
	  t_SubJet_tau4->push_back(SubJet_tau4[i]);   
	}

	t_PV_ndof = PV_ndof;
	t_PV_x = PV_x;
	t_PV_y = PV_y;
	t_PV_z = PV_z;
	t_PV_npvs = PV_npvs;
	t_PV_npvs = PV_npvsGood;
      }
      
      
      tree->Fill();
      clearTreeVectors();


   }
}
