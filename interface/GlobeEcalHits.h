#ifndef GLOBEECALHITS_H
#define GLOBEECALHITS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeEcalHits {
 public:
  
  GlobeEcalHits(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeEcalHits() {};

  void defineBranch(GlobeAnalyzer* ana);
  //bool analyze(const edm::Event&, const edm::EventSetup&, GlobeLeptons*, GlobeElectrons*, GlobeMuons*, GlobePhotons*);
  bool analyze(const edm::Event&, const edm::EventSetup&, GlobeElectrons*, GlobeMuons*, GlobePhotons*);


  // variables

  TClonesArray *ecalhit_p4;
  Int_t ecalhit_n;
  Short_t ecalhit_type[MAX_ECALRECHITS];
  Short_t ecalhit_flag[MAX_ECALRECHITS];
  Float_t ecalhit_time[MAX_ECALRECHITS];
  UInt_t ecalhit_detid[MAX_ECALRECHITS];
  Short_t ecalhit_zside[MAX_ECALRECHITS];
  Short_t ecalhit_ieta[MAX_ECALRECHITS];
  Short_t ecalhit_iphi[MAX_ECALRECHITS];
  Short_t ecalhit_ix[MAX_ECALRECHITS];
  Short_t ecalhit_iy[MAX_ECALRECHITS];
  
 private:
  const char* nome;
  bool doPreshowerHits;
  GlobeCuts *gCUT;
  edm::InputTag ecalHitEBColl, ecalHitEEColl, ecalHitESColl;
  int debug_level;
};


#endif
