#ifndef GLOBEVERTEX_H
#define GLOBEVERTEX_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include <iostream>

class GlobeVertex {
 public:
  
  GlobeVertex(const edm::ParameterSet&, const char*);
  virtual ~GlobeVertex() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t vtx_n;

  TClonesArray *vtx_xyz;
  TClonesArray *vtx_dxdydz;
  TClonesArray *vtx_vectorp3;
  Float_t vtx_scalarpt[MAX_VERTICES];
  Float_t vtx_x2dof[MAX_VERTICES];
  Float_t vtx_ndof[MAX_VERTICES];

  Int_t vtx_ntks[MAX_VERTICES];
  Int_t vtx_tkind[MAX_VERTICES][MAX_VERTEX_TRACKS];
  Float_t vtx_tkweight[MAX_VERTICES][MAX_VERTEX_TRACKS];

 private:
  const char* nome;
  GlobeCuts *gCUT;

  edm::InputTag vertexColl;
  edm::InputTag trackColl;
  int debug_level;
};

#endif
