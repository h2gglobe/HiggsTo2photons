#ifndef GLOBERHO_H
#define GLOBERHO_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "TTree.h"

#include <iostream>

class GlobeRho {
 public:
  
  GlobeRho(const edm::ParameterSet&);
  virtual ~GlobeRho() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  float rho; 

 private:
  edm::InputTag rhoCollection;
  int debug_level;
};

#endif
