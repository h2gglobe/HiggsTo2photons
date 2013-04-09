#ifndef GLOBEGENPARTICLES_H
#define GLOBEGENPARTICLES_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeGenParticles {
 public:
  
  GlobeGenParticles(const edm::ParameterSet&);
  virtual ~GlobeGenParticles() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  //int mother(const HepMC::GenEvent* g, HepMC::GenParticle *p);

  // variables
  Int_t gp_n;

  Short_t gp_status[MAX_GENERATOR];
  Short_t gp_pdgid[MAX_GENERATOR];
  Short_t gp_mother[MAX_GENERATOR];
  
  TClonesArray *gp_p4;
  TClonesArray *gp_vtx;

 private:
  GlobeCuts* gCUT;
  edm::InputTag genParticlesColl;
  int debug_level;
};

#endif
