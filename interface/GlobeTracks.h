#ifndef GLOBETRACKS_H
#define GLOBETRACKS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTrackingParticles.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

// Track Association Methods
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"


#include <iostream>

class GlobeAnalyzer;
class GlobeTrackingParticles;

class GlobeTracks {
 public:
  
  GlobeTracks(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeTracks() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  void GetAssociatedTrackingParticleIndex(const edm::Event&, const edm::EventSetup&, GlobeTrackingParticles*);

  // variables

  Int_t tk_n;
  Int_t tk_nhits[MAX_TRACKS];
  Int_t tk_charge[MAX_TRACKS];
  Int_t tk_nlosthit[MAX_TRACKS];
  Int_t tk_cmsind[MAX_TRACKS];
  Int_t tk_tpind[MAX_TRACKS];
  Float_t tk_chi2[MAX_TRACKS];
  Float_t tk_dof[MAX_TRACKS];
  Float_t tk_d0[MAX_TRACKS];
  Float_t tk_dz[MAX_TRACKS];
  Float_t tk_qoverperr[MAX_TRACKS];
  Float_t tk_pterr[MAX_TRACKS];
  Float_t tk_etaerr[MAX_TRACKS];
  Float_t tk_phierr[MAX_TRACKS];
  Float_t tk_d0err[MAX_TRACKS];
  Float_t tk_dzerr[MAX_TRACKS];
  Int_t tk_hp_nvalid[MAX_TRACKS];
  Int_t tk_hp_nlost[MAX_TRACKS];
  Int_t tk_hp_nvalidpix[MAX_TRACKS];
  Int_t tk_hp_expin[MAX_TRACKS];
  Int_t tk_hp_expout[MAX_TRACKS];
  Int_t tk_algo[MAX_TRACKS];
  Int_t tk_quality[MAX_TRACKS];

  TClonesArray *tk_p4;
  TClonesArray *tk_vtx_pos;

 private:

  std::pair<unsigned int, float> sharedHits(const reco::Track&, const reco::Track&);
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag trackColl;
  edm::InputTag trackColl2;
  edm::InputTag tpColl;
  std::string assocLabel;
  int debug_level;
  bool doAodSim;
};

#endif
