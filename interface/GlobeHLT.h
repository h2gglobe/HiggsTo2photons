#ifndef GLOBEHLT_H
#define GLOBEHLT_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <iostream>
#include <stdint.h>

#include "TClonesArray.h"
#include "TTree.h"

class GlobeHLT {
 public:
  
   GlobeHLT(const edm::ParameterSet&, const char* n="unused");
   virtual ~GlobeHLT() {};

   void defineBranch(TTree* tree);
   bool analyze(const edm::Event&, const edm::EventSetup&);

   std::vector<unsigned short>* hlt_bit;

   Int_t hlt_n;
   std::vector<std::vector<unsigned short> >* hlt_candpath;
   std::vector<std::vector<unsigned short> >* hlt_candpath2;
   std::vector<std::string> *hlt_path_names_HLT;
   std::vector<unsigned int> *hlt_prescale;
  
   TClonesArray* hlt_p4;

   std::vector<std::string> *filter_names_HLT1;
   std::vector<unsigned int>* filter_pass;
   int pass_Mass60_isoiso, pass_Mass60_R9R9, pass_Mass60_mix;
   float Mass60_isoiso, Mass60_isoiso_org;
   int pass_Mass70_isoiso, pass_Mass70_R9R9, pass_Mass70_mix;
   float Mass70_isoiso, Mass70_isoiso_org;

   Int_t PhotonRefs0_n;
   Float_t PhotonRefs0_et[8];
   Float_t PhotonRefs0_eta[8];
   Float_t PhotonRefs0_phi[8];
   Int_t ElectronRefs0_n;
   Float_t ElectronRefs0_et[8];
   Float_t ElectronRefs0_eta[8];
   Float_t ElectronRefs0_phi[8];
   Int_t ElectronRefs1_n;
   Float_t ElectronRefs1_et[8];
   Float_t ElectronRefs1_eta[8];
   Float_t ElectronRefs1_phi[8];
   Int_t ElectronRefs00_n;
   Float_t ElectronRefs00_et[8];
   Float_t ElectronRefs00_eta[8];
   Float_t ElectronRefs00_phi[8];

Int_t ElectronRefs2_n;
   Float_t ElectronRefs2_et[8];
   Float_t ElectronRefs2_eta[8];
   Float_t ElectronRefs2_phi[8];
   Int_t ElectronRefs3_n;
   Float_t ElectronRefs3_et[8];
   Float_t ElectronRefs3_eta[8];
   Float_t ElectronRefs3_phi[8];
 Int_t ElectronRefs4_n;
   Float_t ElectronRefs4_et[8];
   Float_t ElectronRefs4_eta[8];
   Float_t ElectronRefs4_phi[8];
   Int_t ElectronRefs5_n;
   Float_t ElectronRefs5_et[8];
   Float_t ElectronRefs5_eta[8];
   Float_t ElectronRefs5_phi[8];
     
   Int_t PhotonRefs1_n;
   Float_t PhotonRefs1_et[8];
   Float_t PhotonRefs1_eta[8];
   Float_t PhotonRefs1_phi[8];
   Int_t PhotonRefs3_n;
   Float_t PhotonRefs3_et[8];
   Float_t PhotonRefs3_eta[8];
   Float_t PhotonRefs3_phi[8];
   Int_t PhotonRefs4_n;
   Float_t PhotonRefs4_et[8];
   Float_t PhotonRefs4_eta[8];
   Float_t PhotonRefs4_phi[8];
   Int_t PhotonRefs5_n;
   Float_t PhotonRefs5_et[8];
   Float_t PhotonRefs5_eta[8];
   Float_t PhotonRefs5_phi[8];
   Int_t PhotonRefs6_n;
   Float_t PhotonRefs6_et[8];
   Float_t PhotonRefs6_eta[8];
   Float_t PhotonRefs6_phi[8];
   Int_t PhotonRefs8_n;
   Float_t PhotonRefs8_et[8];
   Float_t PhotonRefs8_eta[8];
   Float_t PhotonRefs8_phi[8];
   Int_t PhotonRefs9_n;
   Float_t PhotonRefs9_et[8];
   Float_t PhotonRefs9_eta[8];
   Float_t PhotonRefs9_phi[8];

   Int_t PhotonRefs10_n;
   Float_t PhotonRefs10_et[8];
   Float_t PhotonRefs10_eta[8];
   Float_t PhotonRefs10_phi[8];
   Int_t PhotonRefs11_n;
   Float_t PhotonRefs11_et[8];
   Float_t PhotonRefs11_eta[8];
   Float_t PhotonRefs11_phi[8];





 private:   
   const char* nome;
   GlobeCuts *gCUT;
   edm::InputTag inputTag_;
   edm::InputTag hltTag_;
   int debug_level;

   HLTConfigProvider configProvider;

   std::vector<edm::InputTag> theHLTLabels;
   //bool secondaryTriggerON;
};

#endif
