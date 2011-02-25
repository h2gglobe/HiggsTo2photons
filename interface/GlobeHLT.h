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
  
  std::vector<unsigned short>* hlt1_bit;
  std::vector<unsigned short>* hlt2_bit;

  // MAX_HLT
  Int_t hlt_n;
  std::vector<std::vector<unsigned short> >* hlt_candpath;
  std::vector<std::string> *hlt_path_names_HLT1;
  std::vector<std::string> *hlt_path_names_HLT2;
  
  TClonesArray* hlt_p4;
  
private:   
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag inputTag_;
  edm::InputTag hlt1Tag_;
  edm::InputTag hlt2Tag_;
  int debug_level;
  
  HLTConfigProvider configProvider;
    
  std::vector<edm::InputTag> theHLTLabels;
  bool secondaryTriggerON;
};

#endif
