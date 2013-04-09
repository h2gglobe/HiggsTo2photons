#ifndef GLOBEPILEUP_H
#define GLOBEPILEUP_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TH1D.h"
#include <iostream>
#include <vector>

class GlobePileup {
 public:
  
  GlobePileup(const edm::ParameterSet&);
  virtual ~GlobePileup() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  TH1D* getHisto() { return h1; };
  TH1D* getHistoTrue() { return h2; };

  int pu_n;
  int pu_bunchcrossing;
  float pu_n_true;
  std::vector<float>* pu_zpos;
  std::vector<float>* pu_sumpt_lowpt;
  std::vector<float>* pu_sumpt_highpt;
  std::vector<int>* pu_ntrks_lowpt;
  std::vector<int>* pu_ntrks_highpt;

 private:
   edm::InputTag pileupCollection;
   int debug_level;
   TH1D* h1, *h2;
};

#endif
