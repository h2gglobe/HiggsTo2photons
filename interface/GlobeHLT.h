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
  
  // fullHLT
  Int_t hlt_el_n;
  UInt_t hlt_el_candpath[MAX_HLT];
  Int_t hlt_el_offlineind[MAX_HLT];
  TClonesArray *hlt_el_p4, *hlt_el_tkp4;
  
  Int_t hlt_mu_n;
  UInt_t hlt_mu_candpath[MAX_HLT];
  Int_t hlt_mu_offlineind[MAX_HLT];
  TClonesArray *hlt_mu_p4;
  
  Int_t hlt_ph_n;
  UInt_t hlt_ph_candpath[MAX_HLT];
  Int_t hlt_ph_offlineind[MAX_HLT];
  TClonesArray *hlt_ph_p4;
  
  // new
  // The HLT paths are split into 4 x UInt_t (4 x 32 bit)
  UInt_t hlt1_bit_1;
  UInt_t hlt1_bit_2;
  UInt_t hlt1_bit_3;
  UInt_t hlt1_bit_4;
  UInt_t hlt2_bit_1;
  UInt_t hlt2_bit_2;
  UInt_t hlt2_bit_3;
  UInt_t hlt2_bit_4;
  // MAX_HLT
  Int_t hlt_n;
  // The HLT candpath are split into 4 x UInt_t (4 x 32 bit)
  UInt_t hlt_candpath_1[   MAX_HLT];
  UInt_t hlt_candpath_2[   MAX_HLT];
  UInt_t hlt_candpath_3[   MAX_HLT];
  UInt_t hlt_candpath_4[   MAX_HLT];
  Int_t hlt_id[            MAX_HLT];
  Int_t hlt_jet_offlineind[MAX_HLT];
  
  TClonesArray *hlt_p4;
  
private:   
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag inputTag_;
  edm::InputTag hlt1Tag_;
  edm::InputTag hlt2Tag_;
  edm::InputTag electronColl;
  edm::InputTag muonColl;
  edm::InputTag photonCollStd;
  edm::InputTag jetColl;
  int debug_level;
  bool fullHLT;
  
  // fullHLT
  std::vector<edm::InputTag> theElHLTLabels, theMuHLTLabels, thePhHLTLabels;
  
  // new
  HLTConfigProvider configProvider;
  // The HLT candpaths ans paths are split into 2 x ULong64_t
  //  the corresponding labels are split too
  std::vector<std::string> *hlt_label_names_1;
  std::vector<std::string> *hlt_label_names_2;
  std::vector<std::string> *hlt_label_names_3;
  std::vector<std::string> *hlt_label_names_4;
  std::vector<std::string> *hlt_path_names_HLT1_1;
  std::vector<std::string> *hlt_path_names_HLT1_2;
  std::vector<std::string> *hlt_path_names_HLT1_3;
  std::vector<std::string> *hlt_path_names_HLT1_4;
  std::vector<std::string> *hlt_path_names_HLT2_1;
  std::vector<std::string> *hlt_path_names_HLT2_2;
  std::vector<std::string> *hlt_path_names_HLT2_3;
  std::vector<std::string> *hlt_path_names_HLT2_4;
  
  std::vector<edm::InputTag> theHLTLabels;
  bool secondaryTriggerON;
};

#endif
