#ifndef GLOBEGSFTRACKS_H
#define GLOBEGSFTRACKS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTrackingParticles.h"

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

// Track Association Methods
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"

#include <iostream>

class GlobeTrackingParticles;
class GlobeAnalyzer;

class GlobeGsfTracks {
 public:
  
  GlobeGsfTracks(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeGsfTracks() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  void GetAssociatedTrackingParticleIndex(const edm::Event&, const edm::EventSetup&, GlobeTrackingParticles*);
  
  std::pair<reco::TrackRef,float> getCtfTrackRef(const reco::GsfTrackRef&, edm::Handle<reco::TrackCollection>);
  // variables

  Int_t gsf_tk_n;
  Int_t gsf_tk_nhits[MAX_TRACKS];
  Int_t gsf_tk_charge[MAX_TRACKS];
  Int_t gsf_tk_nlosthit[MAX_TRACKS];
  Int_t gsf_tk_cmsind[MAX_TRACKS];
  Int_t gsf_tk_tpind[MAX_TRACKS];
  Float_t gsf_tk_chi2[MAX_TRACKS];
  Float_t gsf_tk_dof[MAX_TRACKS];
  Float_t gsf_tk_d0[MAX_TRACKS];
  Float_t gsf_tk_dz[MAX_TRACKS];
  Float_t gsf_tk_qoverpouterr[MAX_TRACKS];
  Float_t gsf_tk_qoverpinerr[MAX_TRACKS];
  Float_t gsf_tk_pterr[MAX_TRACKS];
  Float_t gsf_tk_etaerr[MAX_TRACKS];
  Float_t gsf_tk_phierr[MAX_TRACKS];
  Float_t gsf_tk_d0err[MAX_TRACKS];
  Float_t gsf_tk_dzerr[MAX_TRACKS];
  Int_t gsf_tk_hp_nvalid[MAX_TRACKS];
  Int_t gsf_tk_hp_nlost[MAX_TRACKS];
  Int_t gsf_tk_hp_nvalidpix[MAX_TRACKS];
  Int_t gsf_tk_hp_expin[MAX_TRACKS];
  Int_t gsf_tk_hp_expout[MAX_TRACKS];
  Float_t gsf_tk_pin[MAX_TRACKS];
  Float_t gsf_tk_pout[MAX_TRACKS];
  Float_t gsf_tk_fbrem[MAX_TRACKS];
  Int_t gsf_tk_tkind[MAX_TRACKS];
  Float_t gsf_tk_shared[MAX_TRACKS];
  TClonesArray *gsf_tk_p4;
  TClonesArray *gsf_tk_vtx_pos;
  TClonesArray *gsf_tk_poutmode;
  TClonesArray *gsf_tk_pinmode;

 private:
  const char* nome; 
  GlobeCuts *gCUT;
  edm::InputTag trackColl;
  edm::InputTag trackColl2;
  edm::InputTag tpColl;
  edm::InputTag electronColl;
  std::string assocLabel;
  int debug_level;
  bool doAodSim;
  bool storeGsfTracksOnlyIfElectrons;

  edm::ESHandle<MagneticField>                theMagField;
  edm::ESHandle<TrackerGeometry>              trackerHandle_;
  unsigned long long cacheIDGeom_;
  unsigned long long cacheIDTopo_;
  unsigned long long cacheIDTDGeom_;
  unsigned long long cacheIDMagField_;
};

#endif
