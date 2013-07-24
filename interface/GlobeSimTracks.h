#ifndef GLOBESIMTRACKS_H
#define GLOBESIMTRACKS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeSimTracks {
 public:
  
  GlobeSimTracks(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeSimTracks() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables

  Int_t simtrk_n;
  TClonesArray *simtrk_vtx;
  TClonesArray *simtrk_p4;
  Int_t simtrk_pdgid[MAX_SIMTRACKS];
  Int_t simtrk_trkid[MAX_SIMTRACKS];
  Int_t simtrk_mothertrkid[MAX_SIMTRACKS];

  TClonesArray *simvtx;

 private:

  void fillGeantMap();
  SimTrack getNextSimTrack(SimTrack &inputSimtrack);
  SimTrack getSiblingSimTrack(SimTrack &inputSimtrack);
  SimTrack getMotherSimTrack(SimTrack &inputSimtrack);

  const char* nome;
  int debug_level;
  bool track_and_vertex;
  GlobeCuts *gCUT;
  edm::InputTag simtrackColl;
  edm::InputTag simvertexColl; 
  std::map<unsigned, unsigned> geantToIndex_;
  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertices;
};

#endif
