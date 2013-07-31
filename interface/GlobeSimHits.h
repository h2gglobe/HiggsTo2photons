#ifndef GLOBESIMHITS_H
#define GLOBESIMHITS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeSimHits {
 public:
  
  GlobeSimHits(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeSimHits() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables

  Int_t simhit_n;
  TClonesArray *simhit_xyz;
  Float_t simhit_pabs[MAX_SIMHITS];
  Float_t simhit_eloss[MAX_SIMHITS];
  Int_t simhit_subdet[MAX_SIMHITS];
  Int_t simhit_pdgid[MAX_SIMHITS];
  Int_t simhit_trkid[MAX_SIMHITS];
  Int_t simhit_simtrkind[MAX_SIMHITS];

 private:
  const char* nome;
  int debug_level;
  GlobeCuts *gCUT;

  edm::InputTag simhitPixBarrelLowColl;
  edm::InputTag simhitPixEndcapLowColl;
  edm::InputTag simhitTECLowColl;
  edm::InputTag simhitTIBLowColl;
  edm::InputTag simhitTIDLowColl;
  edm::InputTag simhitTOBLowColl;
  edm::InputTag simtrackColl;

  std::vector<PSimHit> theSimHits;
  std::vector<SimTrack> theSimTracks;
};

#endif
