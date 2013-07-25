#ifndef GLOBEGENERATOR_H
#define GLOBEGENERATOR_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeAnalyzer;

class GlobeGenerator {
 public:
  
  GlobeGenerator(const edm::ParameterSet&);
  virtual ~GlobeGenerator() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  int mother(const HepMC::GenEvent* g, HepMC::GenParticle *p);

  // variables
  Int_t gen_n;

  Int_t gen_status[MAX_GENERATOR];
  Int_t gen_pdgid[MAX_GENERATOR];
  Int_t gen_mother[MAX_GENERATOR];
  
  TClonesArray *gen_p4;

 private:
  double etCut_; 
  edm::InputTag generatorColl;
  int debug_level;
  const HepMC::GenEvent* myGenEvent;

};

#endif
