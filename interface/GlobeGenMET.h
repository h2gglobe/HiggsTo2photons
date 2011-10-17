#ifndef GLOBEGENMET_H
#define GLOBEGENMET_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include <iostream>

class GlobeGenMET {
 public:
  
  GlobeGenMET(const edm::ParameterSet&, const char* );
  virtual ~GlobeGenMET() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Float_t met_met;
  Float_t met_phi;
  
 private:
  const char* nome;
  edm::InputTag METColl;
  int debug_level;
};

#endif
