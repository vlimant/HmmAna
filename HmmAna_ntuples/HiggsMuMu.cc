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

      h_mu1pt->Fill((*t_Mu_pt)[0]);
      h_mu2pt->Fill((*t_Mu_pt)[1]);
      h_mu1eta->Fill((*t_Mu_eta)[0]);
      h_mu2eta->Fill((*t_Mu_eta)[1]);
      h_mu1phi->Fill((*t_Mu_phi)[0]);
      h_mu2phi->Fill((*t_Mu_phi)[1]);
      // if (Cut(ientry) < 0) continue;
      if(t_diMuon_mass<120. || t_diMuon_mass>130.){
	h_diMuon_pt->Fill(t_diMuon_pt);
	h_diMuon_eta->Fill(t_diMuon_eta);
	h_diMuon_phi->Fill(t_diMuon_phi);
	h_diMuon_mass->Fill(t_diMuon_mass);
      }
   }
}
