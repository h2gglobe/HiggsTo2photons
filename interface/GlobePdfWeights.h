#ifndef GLOBEPDFWEIGHTS_H
#define GLOBEPDFWEIGHTS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "TTree.h"

class GlobeAnalyzer;

class GlobePdfWeights {
 public:
  
  GlobePdfWeights(const edm::ParameterSet&);
  virtual ~GlobePdfWeights() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  Int_t weight_n;
  Float_t pdf_weights[100];
 private:
  std::vector< edm::InputTag > pdfweightsCollList;
  int debug_level;


};

#endif
