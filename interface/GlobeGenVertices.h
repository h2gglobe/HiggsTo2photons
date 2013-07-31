#ifndef GLOBEGENVERTICES_H
#define GLOBEGENVERTICES_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class GlobeAnalyzer;

class GlobeGenVertices {
 public:
  
  GlobeGenVertices(const edm::ParameterSet&);
  virtual ~GlobeGenVertices() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t gv_n;
  
  TClonesArray* gv_pos;
  TClonesArray* gv_p3;
  
  Float_t gv_sumPtHi[MAX_VERTICES];
  Float_t gv_sumPtLo[MAX_VERTICES];
  Short_t gv_nTkHi[MAX_VERTICES];
  Short_t gv_nTkLo[MAX_VERTICES];

 private:
  edm::InputTag genParticlesColl;
  int debug_level;

};

#endif

