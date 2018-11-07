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
  //  hmm.EventLoop(data, isData);
  hmm.Categorization(data, isData, 100, 150);
  //hmm.GenInfo(data, isData);
  return 0;
}

void HiggsMuMu::GenInfo(const char *data,const char *isData)
{  
   

   double muon_mass = 0.1056583745;
   double el_mass = 0.000511;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      int index_mu1 =0, index_mu2=1;
      vector<int> el;
      vector<int> mu;
      vector<int> extralep;
      el.clear();
      mu.clear();
      extralep.clear();
      TLorentzVector dimu, mu1,mu2;
      
      double evt_wt;
      if(*isData=='T'){evt_wt=1.;}
      else evt_wt=t_genWeight;

      
      int n_GenPart = t_GenPart_pdgId->size();
      for( int i=0; i<n_GenPart; i++){
        for( int j=i+1; j<n_GenPart; j++){
         if(abs((*t_GenPart_pdgId)[i]) == 13 && (*t_GenPart_status)[i] ==1 && abs((*t_GenPart_pdgId)[j]) == 13 && (*t_GenPart_status)[j] ==1){
             mu1.SetPtEtaPhiM((*t_Mu_pt)[index_mu1],(*t_Mu_eta)[index_mu1],(*t_Mu_phi)[index_mu1],muon_mass);
             mu2.SetPtEtaPhiM((*t_Mu_pt)[index_mu2],(*t_Mu_eta)[index_mu2],(*t_Mu_phi)[index_mu2],muon_mass);
             double diMuon_mass = (mu1 + mu2).M();
             if(fabs(diMuon_mass-125.0)<0.5){index_mu1 = i; index_mu2 = j;} 
         }
        }
      }
      /*
      h_gen_mu1mu2dR->Fill(dR,evt_wt);
      h_gen_mu1mu2dPhi->Fill(dPhi,evt_wt);
      h_gen_diMuon_pt->Fill(diMuon_pt,evt_wt);
      h_gen_diMuon_eta->Fill(diMuon_eta,evt_wt);
      */
      for( int i=0; i<n_GenPart; i++){
        if(i!=index_mu1 && i!=index_mu2 && (*t_GenPart_pt)[i] > 10. &&( abs((*t_GenPart_pdgId)[i])== 13 || abs((*t_GenPart_pdgId)[i]) == 11) && (*t_GenPart_status)[i] ==1 ){
        bool overlap=false;
        for(int j=0; j<extralep.size();j++){
           double dRlepH = DeltaR(t_GenPart_eta->at(i),t_GenPart_phi->at(i),t_GenPart_eta->at(extralep.at(j)), t_GenPart_phi->at(extralep.at(j)));
           if(dRlepH<0.1){
              overlap = true;
              break;
           }
           if(overlap) break;
           else extralep.push_back(i);
        }

       
        h_gen_extralep1_pt->Fill(t_GenPart_pt->at(i),evt_wt);
        h_gen_extralep1_eta->Fill(t_GenPart_eta->at(i),evt_wt);
        double dRlepH = DeltaR(t_GenPart_eta->at(i),t_GenPart_phi->at(i),(mu1+mu2).Eta(), (mu1+mu2).Phi());
        h_gen_dRlepH->Fill(dRlepH,evt_wt);
  
      }
     }
     h_gen_extralep->Fill(extralep.size(),evt_wt);
     h_gen_diMuon_m->Fill((mu1 + mu2).M());
  }
}

void HiggsMuMu::Categorization(const char *data,const char *isData, float mlo, float mhi)
{  if (fChain == 0) return;
   cout <<"mass window: "<<mlo<<" - "<<mhi<<endl;
   //TH1F *catyield = new TH1F("h_category_yield","h_category_yield",10,0,10);
   double muon_mass = 0.1056583745;
   double el_mass = 0.000511;
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 2;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      clearTreeVectors();
      double evt_wt;
      if(*isData=='T'){evt_wt=1.;}
      else evt_wt=t_genWeight;
      int index_mu1 =0, index_mu2=1;
      vector<int> el;
      vector<int> mu;
      el.clear();
      /*    mu.clear();

      if(t_Mu_pt->size() > 2){
         bool Event_sel =  false;
         for(int i=0; i<t_Mu_pt->size(); i++){
            if((*t_Mu_pt)[i]>30.){
              for(int j=i+1; j<t_Mu_pt->size(); j++){
                if((*t_Mu_charge)[i]*(*t_Mu_charge)[j]== -1 && (*t_Mu_pt)[j]>20.){
                    Event_sel =true;
                    index_mu1 = i;
                    index_mu2 = j;
                    break;
                 }
              }
            }
            if(Event_sel) break;
         }
	 }
*/
      if(t_Mu_pt->size() > 2){
         for(int i=0; i<t_Mu_pt->size(); i++){
             if(i!=t_mu1 && i!=t_mu2){
                 if(t_Mu_pt->at(i)>10.) mu.push_back(i);
             }
         }  
      }

      if(t_El_pt->size() > 0){
         for(int i=0; i<t_El_pt->size(); i++){
            if(/*t_El_pfRelIso03_all->at(i)<0.15 && */ t_El_pt->at(i)>5 && fabs(t_El_eta->at(i))<2.5 && t_Electron_mvaFall17Iso_WP90->at(i)){ //remove the MVA cut after including the MVA score??
                el.push_back(i);
            }
         }  
      }

      //if(index_mu1!=0 || index_mu2!=1) cout <<"dimuon pair index: "<<index_mu1<<" "<<index_mu2<<endl;
      
      TLorentzVector dimu, mu1,mu2;
      if(t_mu1>1000 || t_mu2>1000)continue;
      mu1.SetPtEtaPhiM((*t_Mu_pt)[t_mu1],(*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],muon_mass);
      mu2.SetPtEtaPhiM((*t_Mu_pt)[t_mu2],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],muon_mass);
      double diMuon_mass = (mu1 + mu2).M(); 
      double diMuon_pt = (mu1 + mu2).Pt();
      double diMuon_eta = (mu1 + mu2).Eta();
      double diMuon_phi = (mu1 + mu2).Phi();
      double dR= DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2]);
      double dEta = (*t_Mu_eta)[t_mu1] - (*t_Mu_eta)[t_mu2];
      double dPhi= DeltaPhi((*t_Mu_phi)[t_mu2],(*t_Mu_phi)[t_mu1]);
      //SR: 120-130, sideband: 110-120, 130-150
      double lepSF;
      if(diMuon_mass>mlo && diMuon_mass<mhi){
          double binv = catyield->GetBinContent(10);
          binv = binv + t_genWeight;
	  if(*isData=='F'){
	    if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_TRIG)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
	    
	    else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_TRIG)[t_mu2]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
	    
	    
	    evt_wt*=lepSF;
	  }
          catyield->SetBinContent(10,binv);
          h_diMuon_mass_cat->Fill(diMuon_mass,evt_wt);
        
          run  =  t_run;
          event = t_event;
          lumi = t_luminosityBlock;
          genWeight = evt_wt;
          Higgs_mass = diMuon_mass;
          Higgs_pt = diMuon_pt;
          Higgs_eta = diMuon_eta; 
          //ttH
          if(t_nbJet>0){

	    //	    if(t_nJet!=t_nbJet)cout<<evt_wt<<" :";
	    if(*isData=='F'){
	      /*for(int k=0;k<t_nJet;k++){
		//		if(t_nJet!=t_nbJet)cout<<fabs(t_nJet-t_nbJet)<<endl;
		if((*t_Jet_btagDeepB)[k]>0.4941) evt_wt*=(*t_bJet_SF)[k];
		else evt_wt*=(1-(*t_bJet_SF)[k]);
		}*/
	    }
	    //	    if(t_nJet!=t_nbJet)cout<<evt_wt<<endl;
	    if(el.size()> 0  || mu.size()>0){
	      cat_index = 1;
	      double binv = catyield->GetBinContent(1);
	      binv = binv + t_genWeight;
	      catyield->SetBinContent(1,binv);
	      if(*isData=='F')h_diMuon_mass_ttHLep->Fill(diMuon_mass,evt_wt);
	      if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_ttHLep->Fill(diMuon_mass,evt_wt);
	      else if(*isData=='F')h_diMuon_mass_110To140_ttHLep->Fill(diMuon_mass,evt_wt);
	      if(diMuon_mass<120. || diMuon_mass>130.){
		h_mu1mu2dR_ttHLep->Fill(dR,evt_wt);
		h_mu1mu2dPhi_ttHLep->Fill(dPhi,evt_wt);
		h_mu1mu2dEta_ttHLep->Fill(dEta,evt_wt);
		h_diMuon_pt_ttHLep->Fill(diMuon_pt,evt_wt);
		h_diMuon_eta_ttHLep->Fill(diMuon_eta,evt_wt);
		h_diMuon_phi_ttHLep->Fill(diMuon_phi,evt_wt);
		//h_Nbjet_ttHLep->Fill(t_nbJet);
		h_leading_bJet_pt_ttHLep->Fill((*t_bJet_pt)[0],evt_wt);
		h_leading_bJet_eta_ttHLep->Fill((*t_bJet_eta)[0],evt_wt);
		h_leading_bJet_phi_ttHLep->Fill((*t_bJet_phi)[0],evt_wt);
	      }
	    }
              else if(t_nJet>4){
                  cat_index = 2;
                  double binv = catyield->GetBinContent(2);
                  binv = binv + t_genWeight;
                  catyield->SetBinContent(2,binv);
                  if(*isData=='F')h_diMuon_mass_ttHHad->Fill(diMuon_mass,evt_wt);
		  if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_ttHHad->Fill(diMuon_mass,evt_wt);
		  else if(*isData=='F')h_diMuon_mass_110To140_ttHHad->Fill(diMuon_mass,evt_wt);
		  if(diMuon_mass<120. || diMuon_mass>130.){
		    h_mu1mu2dR_ttHHad->Fill(dR,evt_wt);
		    h_mu1mu2dPhi_ttHHad->Fill(dPhi,evt_wt);
		    h_mu1mu2dEta_ttHHad->Fill(dEta,evt_wt);
		    h_diMuon_pt_ttHHad->Fill(diMuon_pt,evt_wt);
		    h_diMuon_eta_ttHHad->Fill(diMuon_eta,evt_wt);
		    h_diMuon_phi_ttHHad->Fill(diMuon_phi,evt_wt);
		    h_leading_bJet_pt_ttHHad->Fill((*t_bJet_pt)[0],evt_wt);
		    h_leading_bJet_eta_ttHHad->Fill((*t_bJet_eta)[0],evt_wt);
		    h_leading_bJet_phi_ttHHad->Fill((*t_bJet_phi)[0],evt_wt);
		  }
              }
              else{
                  cat_index = 3;
                  double binv = catyield->GetBinContent(3);
                  binv = binv + t_genWeight;
                  catyield->SetBinContent(3,binv);
                  if(*isData=='F')h_diMuon_mass_ttHLoose->Fill(diMuon_mass,evt_wt);
		  if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_ttHLoose->Fill(diMuon_mass,evt_wt);
		  else if(*isData=='F')h_diMuon_mass_110To140_ttHLoose->Fill(diMuon_mass,evt_wt);
		  if(diMuon_mass<120. || diMuon_mass>130.){
		    h_mu1mu2dR_ttHLoose->Fill(dR,evt_wt);
		    h_mu1mu2dPhi_ttHLoose->Fill(dPhi,evt_wt);
		    h_mu1mu2dEta_ttHLoose->Fill(dEta,evt_wt);
		    h_diMuon_pt_ttHLoose->Fill(diMuon_pt,evt_wt);
		    h_diMuon_eta_ttHLoose->Fill(diMuon_eta,evt_wt);
		    h_diMuon_phi_ttHLoose->Fill(diMuon_phi,evt_wt);
		    h_leading_bJet_pt_ttHLoose->Fill((*t_bJet_pt)[0],evt_wt);
		    h_leading_bJet_eta_ttHLoose->Fill((*t_bJet_eta)[0],evt_wt);
		    h_leading_bJet_phi_ttHLoose->Fill((*t_bJet_phi)[0],evt_wt);
		  }
              }
              cattree->Fill();
	  }
          //ZH mm
          else if(mu.size()>1){
              bool isOS=false;
              for(int k=0; k<mu.size(); k++){
                   int m1_chg = t_Mu_charge->at(mu.at(k));
                   for(int j=k+1; j<mu.size(); j++){
                       int m2_chg = t_Mu_charge->at(mu.at(j));
                       if(m1_chg*m2_chg==-1){
                              TLorentzVector extra_m1, extra_m2;
                              extra_m1.SetPtEtaPhiM((*t_Mu_pt)[mu.at(k)],(*t_Mu_eta)[mu.at(k)],(*t_Mu_phi)[mu.at(k)],muon_mass);             
                              extra_m2.SetPtEtaPhiM((*t_Mu_pt)[mu.at(j)],(*t_Mu_eta)[mu.at(j)],(*t_Mu_phi)[mu.at(j)],muon_mass);             
                              double dimu_mass = (extra_m1+extra_m2).M();
                              if(fabs(dimu_mass-91.1876) < 20.){
                                    isOS=true;
                                    l1_index->push_back(mu.at(k));
                                    l2_index->push_back(mu.at(j));
                              }
                       }
                   }
              }

              if(isOS){
                  cat_index = 4;
                  double binv = catyield->GetBinContent(4);
                  binv = binv + t_genWeight;
                  catyield->SetBinContent(4,binv);
                  h_diMuon_mass_ZHll->Fill(diMuon_mass,evt_wt);
                  cattree->Fill();
              }
          }
          //ZH ee
          else if(t_El_pt->size()>1){
              bool isOS=false; 
              for(int k=0; k<el.size(); k++){
                   int e1_chg = t_El_charge->at(el.at(k));
                   for(int j=k+1; j<el.size(); j++){
                       int e2_chg = t_El_charge->at(el.at(j));
                       if(e1_chg*e2_chg==-1){
                              TLorentzVector extra_e1, extra_e2;
                              extra_e1.SetPtEtaPhiM((*t_El_pt)[el.at(k)],(*t_El_eta)[el.at(k)],(*t_El_phi)[el.at(k)],el_mass);
                              extra_e2.SetPtEtaPhiM((*t_El_pt)[el.at(j)],(*t_El_eta)[el.at(j)],(*t_El_phi)[el.at(j)],el_mass);
                              double diel_mass = (extra_e1+extra_e2).M();
                              if(fabs(diel_mass-91.1876) < 20.){
                                    isOS=true;
                                    l1_index->push_back(el.at(k));
                                    l2_index->push_back(el.at(j));
                              }
                       }
                   }
              }
         
              if(isOS){
                  cat_index = 4;
                  double binv = catyield->GetBinContent(4);
                  binv = binv + t_genWeight;
                  catyield->SetBinContent(4,binv);
                  h_diMuon_mass_ZHll->Fill(diMuon_mass,evt_wt);
                  cattree->Fill();
              }
          }
          //WH, W->ev
          else if(el.size()>0){

              double binv = catyield->GetBinContent(5);
              binv = binv + t_genWeight;
              catyield->SetBinContent(5,binv);
	      double dRlepH = DeltaR(t_El_eta->at(el.at(0)),t_El_phi->at(el.at(0)),(mu1+mu2).Eta(), (mu1+mu2).Phi());
	
	      if(*isData=='F')h_diMuon_mass_WHlv->Fill(diMuon_mass,evt_wt);
	      if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_WHlv->Fill(diMuon_mass,evt_wt);
	      else if(*isData=='F')h_diMuon_mass_110To140_WHlv->Fill(diMuon_mass,evt_wt);
              //h_mu1mu2dEta->Fill(dEta,evt_wt); 
	      if(diMuon_mass<120. || diMuon_mass>130.){
                h_MET_pt_WHTolv->Fill(t_MET_pt,evt_wt);
		h_mu1mu2dR_WHTolv->Fill(dR,evt_wt);
		h_mu1mu2dPhi_WHTolv->Fill(dPhi,evt_wt);
		h_mu1mu2dEta_WHTolv->Fill(dEta,evt_wt);
		h_diMuon_pt_WHTolv->Fill(diMuon_pt,evt_wt);
		h_diMuon_eta_WHTolv->Fill(diMuon_eta,evt_wt);
		h_diMuon_phi_WHTolv->Fill(diMuon_phi,evt_wt);
		h_extralep1_pt->Fill(t_El_pt->at(el.at(0)),evt_wt);
		h_extralep1_eta->Fill(t_El_eta->at(el.at(0)),evt_wt);
		h_dRlepH->Fill(dRlepH,evt_wt); 
		h_extralep_Electron_mvaFall17Iso->Fill(t_Electron_mvaFall17Iso->at(el.at(0)),evt_wt);
		
	      }
              cat_index = 5;
              /*MET_pt = t_MET_pt;
              MET_phi = t_MET_phi;
              extralep_pfRelIso03 = t_El_pfRelIso03_all->at(el.at(0));
              extralep_pt = t_El_pt->at(el.at(0));
              extralep_eta = t_El_eta->at(el.at(0)); 
	      dRlepHiggs = dRlepH;
              dRmm = dR;
              dEtamm = dEta;
              dPhimm = dPhi;*/
              cattree->Fill();
          }
          //WH, W->mv
          else if(mu.size()>0){
          //else if(t_Mu_pt->size()==3 && t_Mu_pt->at(2)>10. && fabs(t_Mu_eta->at(2))<2.4){
	    if(*isData=='F')evt_wt*=t_Mu_EffSF_ID->at(mu.at(0))*t_Mu_EffSF_ISO->at(mu.at(0));
	    double binv = catyield->GetBinContent(6);
	    binv = binv + t_genWeight;
	    catyield->SetBinContent(6,binv);
	    double dRlepH = DeltaR(t_Mu_eta->at(mu.at(0)),t_Mu_phi->at(mu.at(0)),(mu1+mu2).Eta(), (mu1+mu2).Phi());
             
	    if(*isData=='F')h_diMuon_mass_WHlv->Fill(diMuon_mass,evt_wt);
	    if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_WHlv->Fill(diMuon_mass,evt_wt);
	    else if(*isData=='F')h_diMuon_mass_110To140_WHlv->Fill(diMuon_mass,evt_wt);
	    if(diMuon_mass<120. || diMuon_mass>130.){
              h_MET_pt_WHTolv->Fill(t_MET_pt,evt_wt);
	      h_mu1mu2dR_WHTolv->Fill(dR,evt_wt);
	      h_mu1mu2dPhi_WHTolv->Fill(dPhi,evt_wt);
	      h_mu1mu2dEta_WHTolv->Fill(dEta,evt_wt);
	      h_diMuon_pt_WHTolv->Fill(diMuon_pt,evt_wt);
	      h_diMuon_eta_WHTolv->Fill(diMuon_eta,evt_wt);
	      h_diMuon_phi_WHTolv->Fill(diMuon_phi,evt_wt);
              h_extralep1_pt->Fill(t_Mu_pt->at(mu.at(0)),evt_wt);
              h_extralep1_eta->Fill(t_Mu_eta->at(mu.at(0)),evt_wt);
              h_dRlepH->Fill(dRlepH,evt_wt);
	    }
              //run  =  t_run;
              //event =  t_event;
              //genWeight = evt_wt;
              cat_index = 6;
              /*MET_pt = t_MET_pt;
              MET_phi = t_MET_phi;
              extralep_pt = t_Mu_pt->at(mu.at(0));
              extralep_eta = t_Mu_eta->at(mu.at(0));
              dRlepHiggs = dRlepH;
              dRmm = dR;
              dEtamm = dEta;
              dPhimm = dPhi;*/
              cattree->Fill();
          }
          //VBF
          else if(t_diJet_mass>400.){
              cat_index = 7;
              double binv = catyield->GetBinContent(7);
              binv = binv + t_genWeight;
              catyield->SetBinContent(7,binv);
              if(*isData=='F')h_diMuon_mass_VBF->Fill(diMuon_mass,evt_wt);
	      if(diMuon_mass<120. || diMuon_mass>130.)h_diMuon_mass_110To140_VBF->Fill(diMuon_mass,evt_wt);
	      else if(*isData=='F')h_diMuon_mass_110To140_VBF->Fill(diMuon_mass,evt_wt);
	      if(diMuon_mass<120. || diMuon_mass>130.){
		h_diMuon_pt_VBF->Fill(diMuon_pt,evt_wt);
		h_diMuon_eta_VBF->Fill(diMuon_eta,evt_wt);
		h_diMuon_phi_VBF->Fill(diMuon_phi,evt_wt);
		h_mu1mu2dR_VBF->Fill(dR,evt_wt);
		h_mu1mu2dPhi_VBF->Fill(dPhi,evt_wt);
		h_mu1mu2dEta_VBF->Fill(dEta,evt_wt);

		h_dijet_pt_VBF->Fill(t_diJet_pt,evt_wt);
		h_dijet_eta_VBF->Fill(t_diJet_eta,evt_wt);
		h_dijet_phi_VBF->Fill(t_diJet_phi,evt_wt);
		h_Mjj_VBF->Fill(t_diJet_mass,evt_wt);
		h_dijet_dEta_VBF->Fill((*t_Jet_eta)[0]-(*t_Jet_eta)[1],evt_wt);
	      }
	      dRmm = dR;
              dEtamm = dEta;
              dPhimm = dPhi;
	      M_jj=t_diJet_mass;
	      dEta_jj=(*t_Jet_eta)[0]-(*t_Jet_eta)[1];
	      Zep=diMuon_eta-0.5*((*t_Jet_eta)[0]+(*t_Jet_eta)[1])/fabs((*t_Jet_eta)[0]-(*t_Jet_eta)[1]);
	      
	      double dr[4];
	      dr[0]=DeltaR((*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
	      dr[1] = DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
	      dr[2]=DeltaR((*t_Mu_eta)[t_mu2],(*t_Mu_phi)[t_mu2],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
	      dr[3]=DeltaR((*t_Mu_eta)[t_mu1],(*t_Mu_phi)[t_mu1],(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
	      
	      for (int c = 0 ; c < 3; c++)
		{
		  for (int d = 0 ; d < 3 - c; d++)
		    {
		      if (dr[d] > dr[d+1]) /* For decreasing order use < */
			{
			  double swap = dr[d];
			  dr[d]   = dr[d+1];
			  dr[d+1] = swap;
			}
		    }
		}
	      dRmin_mj=dr[0];
	      dRmax_mj=dr[3];

	      dr[0]=DeltaR(diMuon_eta,diMuon_phi,(*t_Jet_eta)[0],(*t_Jet_phi)[0]);
	      dr[1]=DeltaR(diMuon_eta,diMuon_phi,(*t_Jet_eta)[1],(*t_Jet_phi)[1]);
	      
	      if(dr[0]>dr[1]){
		dRmin_mmj=dr[1];
		dRmax_mmj=dr[0];
	      }
	      else{
		dRmin_mmj=dr[0];
                dRmax_mmj=dr[1];
	      }
	      
	      dPhijj=DeltaPhi((*t_Jet_phi)[1],(*t_Jet_phi)[0]);
	      leadingJet_pt=(*t_Jet_pt)[0];
	      subleadingJet_pt=(*t_Jet_pt)[1];
	      leadingJet_eta = (*t_Jet_eta)[0];
	      subleadingJet_eta = (*t_Jet_eta)[1];

	      cthetaCS=2*(mu2.E()*mu1.Pz()-mu1.E()*mu2.Pz())/(diMuon_mass*sqrt(pow(diMuon_mass,2)+pow(diMuon_pt,2)));
	      cout<<event<<" "<<cthetaCS<<" "<<Zep<<endl;
	      cattree->Fill();
	      if(diMuon_mass<120. || diMuon_mass>130.){
		h_dijet_dPhijj_VBF->Fill(dPhijj,evt_wt);
		h_jet1_pt_VBF->Fill(leadingJet_pt,evt_wt);
		h_jet2_pt_VBF->Fill(subleadingJet_pt,evt_wt);
		h_jet1_eta_VBF->Fill(leadingJet_eta,evt_wt);
		h_jet2_eta_VBF->Fill(subleadingJet_eta,evt_wt);
		h_cthetaCS_VBF->Fill(cthetaCS,evt_wt);
		h_dRmin_mj_VBF->Fill(dRmin_mj,evt_wt);
		h_dRmax_mj_VBF->Fill(dRmax_mj,evt_wt);
		h_dRmin_mmj_VBF->Fill(dRmin_mmj,evt_wt);
                h_dRmax_mmj_VBF->Fill(dRmax_mmj,evt_wt);
		h_Zep_VBF->Fill(Zep,evt_wt);
	      }
	  }
	  
          //VH, had
          else if(t_diJet_mass>60. && t_diJet_mass<110){
              cat_index = 8;
              double binv = catyield->GetBinContent(8);
              binv = binv + t_genWeight;
              catyield->SetBinContent(8,binv);
              h_diMuon_mass_VHHad->Fill(diMuon_mass,evt_wt);
              cattree->Fill();
          }
          //ggH
          else{
              cat_index = 9;
              double binv = catyield->GetBinContent(9);
              binv = binv + t_genWeight;
              catyield->SetBinContent(9,binv);
              h_diMuon_mass_ggH->Fill(diMuon_mass,evt_wt);
              cattree->Fill();
          }
      }
   }
}

void HiggsMuMu::EventLoop(const char *data,const char *isData)
{  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<jentry<<endl;
      double evt_wt;
      if(*isData=='T'){evt_wt=1.;}
      else evt_wt=t_genWeight;
      TLorentzVector m1,m2,m12;
      bool olp=false;
      int mu2idx=1;
      if((*t_Mu_charge)[0]*(*t_Mu_charge)[1]==-1){
	//	m1.SetPtEtaPhiM((*t_Mu_pt)[0],(*t_Mu_eta)[0],(*t_Mu_phi)[0],0.1056583745);
	//m2.SetPtEtaPhiM((*t_Mu_pt)[1],(*t_Mu_eta)[1],(*t_Mu_phi)[1],0.1056583745);
	m1.SetPtEtaPhiM((*t_Mu_pt)[0],(*t_Mu_eta)[0],(*t_Mu_phi)[0],(*t_Mu_mass)[0]);
	m2.SetPtEtaPhiM((*t_Mu_pt)[1],(*t_Mu_eta)[1],(*t_Mu_phi)[1],(*t_Mu_mass)[1]);
	m12=m1+m2;
	olp=true;
      }
      else{
	for(int i=2;i<t_Mu_pt->size();i++){
	  if((*t_Mu_pt)[i]>20. && (*t_Mu_charge)[0]*(*t_Mu_charge)[i]==-1){
	    //	    m1.SetPtEtaPhiM((*t_Mu_pt)[0],(*t_Mu_eta)[0],(*t_Mu_phi)[0],0.1056583745);
	    //m2.SetPtEtaPhiM((*t_Mu_pt)[i],(*t_Mu_eta)[i],(*t_Mu_phi)[i],0.1056583745);
	    m1.SetPtEtaPhiM((*t_Mu_pt)[0],(*t_Mu_eta)[0],(*t_Mu_phi)[0],(*t_Mu_mass)[0]);
	    m2.SetPtEtaPhiM((*t_Mu_pt)[1],(*t_Mu_eta)[1],(*t_Mu_phi)[1],(*t_Mu_mass)[1]);
	    m12=m1+m2;
	    olp=true;
	    mu2idx=i;
	    break;
	  }
	}
      }
      //cout<<olp<<endl;


      double lepSF;
      if(olp){

	if(*isData=='F'){
	  if(t_index_trigm_mu==t_mu1) lepSF = (*t_Mu_EffSF_TRIG)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
	  
	  else if( t_index_trigm_mu==t_mu2) lepSF = (*t_Mu_EffSF_TRIG)[t_mu2]*(*t_Mu_EffSF_ID)[t_mu1]*(*t_Mu_EffSF_ISO)[t_mu1]*(*t_Mu_EffSF_ID)[t_mu2]*(*t_Mu_EffSF_ISO)[t_mu2];
	    
	    
	  evt_wt*=lepSF;
	}
	//cout<<evt_wt<<endl;
	//if(t_diMuon_mass<120. || t_diMuon_mass>130.)h_diMuon_mass->Fill(t_diMuon_mass,evt_wt);
	if(*isData=='T'){
	  if(m12.M()<120. || m12.M()>130.)
	    h_diMuon_mass->Fill(m12.M(),evt_wt);
	}
	else h_diMuon_mass->Fill(m12.M(),evt_wt);
	if(m12.M()>120. && m12.M()<130.)h_diMuon_mass_SR->Fill(m12.M(),evt_wt);
	if((m12.M()>110. && m12.M()<120.) || (m12.M()>130. && m12.M()<150.)){
	  h_mu1pt->Fill((*t_Mu_pt)[0],evt_wt);
	  h_mu2pt->Fill((*t_Mu_pt)[mu2idx],evt_wt);
	  h_mu1eta->Fill((*t_Mu_eta)[0],evt_wt);
	  h_mu2eta->Fill((*t_Mu_eta)[mu2idx],evt_wt);
	  h_mu1phi->Fill((*t_Mu_phi)[0],evt_wt);
	  h_mu2phi->Fill((*t_Mu_phi)[mu2idx],evt_wt);
	  double dR= DeltaR((*t_Mu_eta)[0],(*t_Mu_phi)[0],(*t_Mu_eta)[mu2idx],(*t_Mu_phi)[mu2idx]);
	  h_mu1mu2dR->Fill(dR,evt_wt);
	  double dPhi= DeltaPhi((*t_Mu_phi)[mu2idx],(*t_Mu_phi)[0]);
	  h_mu1mu2dPhi->Fill(dPhi,evt_wt);
	  // if (Cut(ientry) < 0) continue;
	
	  //cout<<"done till dPhi\n";
		
	  h_diMuon_pt->Fill(m12.Pt(),evt_wt);
	  h_diMuon_eta->Fill(m12.Eta(),evt_wt);
	  h_diMuon_phi->Fill(m12.Phi(),evt_wt);
	  h_diMuon_mass_110To150->Fill(m12.M(),evt_wt);
	  if(m12.M()>110.&& m12.M()<120.)h_diMuon_mass_110To120->Fill(m12.M(),evt_wt);
	  if(m12.M()>130. && m12.M()<150.)h_diMuon_mass_130To150->Fill(m12.M(),evt_wt);
	  //cout<<"done all dimu\n";	
	  h_MET_pt->Fill(t_MET_pt,evt_wt);
	  h_METphi->Fill(t_MET_phi,evt_wt);
	  h_MET_sumEt->Fill(t_MET_sumEt,evt_wt);

	  h_Njet->Fill(t_nJet,evt_wt);
	  h_Nbjet->Fill(t_nbJet,evt_wt);
	  //cout<<t_nJet<<endl;
	  if(t_nJet>0){
	    //cout<<(*t_Jet_pt)[0]<<endl;
	    h_j1pt->Fill((*t_Jet_pt)[0],evt_wt);
	    //cout<<"pt\n";
	    h_j1eta->Fill((*t_Jet_eta)[0],evt_wt);
	    h_j1phi->Fill((*t_Jet_phi)[0],evt_wt);
	  }
	  //cout<<"done with jet 1\n";
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
      //cout<<"ok done\n";
   }
}
