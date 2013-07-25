#ifndef GLOBEGENJETS_H
#define GLOBEGENJETS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeAnalyzer;

class GlobeGenJets {
 public:
  
  GlobeGenJets(const edm::ParameterSet&, const char*);
  virtual ~GlobeGenJets() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t genjet_n;

  Float_t genjet_em[MAX_GENJETS]; 
  Float_t genjet_had[MAX_GENJETS]; 
  Float_t genjet_inv[MAX_GENJETS]; 
  Float_t genjet_aux[MAX_GENJETS];
  
  TClonesArray *genjet_p4;

 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag jetColl;
  int debug_level;
};

#endif
