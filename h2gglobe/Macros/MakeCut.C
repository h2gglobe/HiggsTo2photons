#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

void makeNewTree(TTree *oldTree, TTree *newTree, TFile *outFile){
  
  float bdtoutput, deltaMOverM, wt;
  oldTree->SetBranchAddress("bdtoutput",&bdtoutput);
  oldTree->SetBranchAddress("deltaMOverM",&deltaMOverM);
  oldTree->SetBranchAddress("wt",&wt);
  newTree->Branch("bdtoutput",&bdtoutput,"bdtoutput/F");
  newTree->Branch("deltaMOverM",&deltaMOverM,"deltaMOverM/F");
  newTree->Branch("wt",&wt,"wt/F");

  string name = string(newTree->GetName());
  TH1F *h_bdt = new TH1F(Form("h_bdt_%s",name.c_str()),Form("h_bdt_%s",name.c_str()),100,-1,1.);
  TH1F *h_dM = new TH1F(Form("h_dm_%s",name.c_str()),Form("h_dm_%s",name.c_str()),100,-0.1,0.1);
  TH2F *scat = new TH2F(Form("scat_%s",name.c_str()),Form("scat_%s",name.c_str()),100,-1.,1.,100,-.1,0.1);
  TH2F *scat_abs = new TH2F(Form("scat_abs_%s",name.c_str()),Form("scat_abs_%s",name.c_str()),100,-1.,1.,100,0,0.1);
  
  for (int i=0; i<oldTree->GetEntries(); i++){
    oldTree->GetEntry(i);
    if (bdtoutput > -0.5 && deltaMOverM <= 0.02) {
      newTree->Fill();
      h_bdt->Fill(bdtoutput);
      h_dM->Fill(deltaMOverM);
      scat->Fill(bdtoutput,deltaMOverM);
      scat_abs->Fill(bdtoutput,fabs(deltaMOverM));
    }
  }
  outFile->cd();
  h_bdt->Write();
  h_dM->Write();
  scat->Write();
  scat_abs->Write();

}

void MakeCut(){
  
  TFile *inFile = TFile::Open("../../PhotonAnalysis_scripts/TMVA_input_CMS-HGG_mit_training_07_01_12.root");
  TFile *outFile = new TFile("TMVA_input_CMS-HGG_mit_training_07_01_12_slim.root","RECREATE");

  const int nTrees=4;
  string oldTreeNames[nTrees] = {"sig_train_123.0","sig_test_123.0","bkg_train_2pt_123.0","bkg_test_2pt_123.0"};
  string newTreeNames[nTrees] = {"sig_train","sig_test","bkg_train","bkg_test"};

  for (int n=0; n<nTrees; n++){
    TTree *oldTree = (TTree*)inFile->Get(oldTreeNames[n].c_str());
    TTree *newTree = new TTree(newTreeNames[n].c_str(),newTreeNames[n].c_str());
    makeNewTree(oldTree,newTree,outFile);
    outFile->cd();
    newTree->Write();
  }


}
