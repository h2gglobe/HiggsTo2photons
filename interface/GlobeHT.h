#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEHT_H
#define GLOBEHT_H

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMET.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCaloTowers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeHT {
 public:
  
  GlobeHT(const edm::ParameterSet&);
  virtual ~GlobeHT() {};

  void defineBranch(GlobeAnalyzer* ana);
  void fillCaloTowerHT(GlobeMET * theMET, GlobeCaloTowers * theCaloTowers);
  void fillTrackHT(const edm::Event&);
  void fillLeptonHT(GlobeJets * theJets, GlobeMET * theMET, GlobeLeptons * theLeptons);

  //variables
  Float_t ht_nomet25;
  Float_t ht_nomet35;
  Float_t ht_nomet50;
  Float_t ht_trk;
  Float_t ht_25;
  Float_t ht_35;
  Float_t ht_50;
  TVector3 *ht_trkvec;

  Int_t ht_2lpt_n;
  Int_t ht_2lpt_inds[MAX_HT2][2];
  Float_t ht_2lpt25[MAX_HT2];
  Float_t ht_2lpt35[MAX_HT2];
  Float_t ht_2lpt50[MAX_HT2];

  Int_t ht_4lpt_n;
  Int_t ht_4lpt_inds[MAX_HT4][4];
  Float_t ht_4lpt25[MAX_HT4];
  Float_t ht_4lpt35[MAX_HT4];
  Float_t ht_4lpt50[MAX_HT4];


 private:
  int debug_level;
  edm::InputTag trackColl;
};

#endif
