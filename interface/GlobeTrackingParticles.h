#ifndef GLOBETRACKINGPARTICLES_H
#define GLOBETRACKINGPARTICLES_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// Tracking Object Definitions
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Magnetic Field and track propagation 
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Track Association Methods
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include <iostream>

using namespace HepMC;
using namespace std;

class GlobeTracks;
class GlobeAnalyzer;

class GlobeTrackingParticles {
 public:

  typedef edm::RefVector<edm::HepMCProduct, HepMC::GenParticle > GenParticleRefVector;
  typedef edm::Ref<edm::HepMCProduct, HepMC::GenParticle >       GenParticleRef;
  
  GlobeTrackingParticles(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeTrackingParticles() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&, GlobeTracks *);

  // variables
  Int_t tp_n;
  TClonesArray* tp_vtx;
  TClonesArray* tp_p4;
  TClonesArray* tv_xyz;
  Int_t tp_pdgid[MAX_TRACKINGPARTICLES];
  Int_t tp_motherid[MAX_TRACKINGPARTICLES];
  Int_t tp_charge[MAX_TRACKINGPARTICLES];
  Int_t tp_tkind[MAX_TRACKINGPARTICLES];
  Int_t tp_genind[MAX_TRACKINGPARTICLES];
  Int_t tp_cmsind[MAX_TRACKINGPARTICLES]; //not filled
  Double_t tp_d0[MAX_TRACKINGPARTICLES];
  Double_t tp_dz[MAX_TRACKINGPARTICLES];
 
 protected:
  
  void CalcTrackingParticleImpactParam(double &d0, double &dz, TrackingParticleRef, const MagneticField *bf);
  Int_t GetAssociatedRecoTrackIndex(TrackingParticleRef tp, const edm::Handle<edm::View<reco::Track> >& tkH, const reco::SimToRecoCollection &simRecColl, GlobeTracks*);
  Int_t GetAssociatedGenParticleIndex(TrackingParticleRef tp, const edm::Handle<edm::HepMCProduct>& genEvtH);
  Int_t GetMotherTrackingParticle(TrackingParticleRef tp) {return GetMotherTrackingParticle(*tp);}
  Int_t GetMotherTrackingParticle(const TrackingParticle&);

 private:
    
  const char* nome;
  int debug_level;
  edm::InputTag tpColl;
  edm::InputTag tvColl;
  edm::InputTag trackColl;
  std::string assocLabel;
  GlobeCuts *gCUT;

  // need this stuff for mother (hopefully only temporary)
  edm::InputTag simtrackColl;
  edm::InputTag simvertexColl;
  edm::InputTag generatorColl;
  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertices;

  Int_t GetMotherSimTrack(const SimTrack&);
  Int_t GetMotherGenParticle(const HepMC::GenParticle&);
  
};

#endif
