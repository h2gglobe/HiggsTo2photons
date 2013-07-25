#ifndef GLOBEHCAL_H
#define GLOBEHCAL_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeHcal {
 public:
  
  GlobeHcal(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeHcal() {};

  void defineBranch(GlobeAnalyzer* ana);
  //bool analyze(const edm::Event&, const edm::EventSetup&,GlobeLeptons*,GlobeElectrons*,GlobeMuons*,GlobePhotons*,GlobeTracks*);
  bool analyze(const edm::Event&, const edm::EventSetup&, GlobeElectrons*, GlobeMuons*, GlobePhotons*);

  // variables

  TClonesArray *hc_p4;
  Int_t hc_type[MAX_HCALHITS];
  Int_t hc_n;

 private:
  const char* nome;
  bool doHFHcal;
  GlobeCuts *gCUT;
  edm::InputTag hcalBEColl;
  edm::InputTag hcalFColl;
  edm::InputTag hcalHoColl; 
  int debug_level;
};


#endif
