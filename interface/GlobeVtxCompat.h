#ifndef GLOBEVTXCOMPAT_H
#define GLOBEVTXCOMPAT_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>
#include <vector>

struct lepton { int trkind; int lepind; };

class GlobeVtxCompat {
 public:
  
  GlobeVtxCompat(const edm::ParameterSet&);
  virtual ~GlobeVtxCompat() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&, GlobeLeptons*, GlobeElectrons*, GlobeMuons*, GlobeTracks*);

  // variables
  Int_t vtxcomp_n;
  Int_t vtxcomp_lepn[MAX_VTXCOMP];
  Int_t vtxcomp_pdgid[MAX_VTXCOMP][4];
  Int_t vtxcomp_lepind[MAX_VTXCOMP][4];

  Bool_t  vtxcomp_klmn_isvalid[MAX_VTXCOMP];  
  Float_t vtxcomp_klmn_x2dof[MAX_VTXCOMP];    
  Int_t   vtxcomp_klmn_ndof[MAX_VTXCOMP];     
  Float_t vtxcomp_klmn_x2prob[MAX_VTXCOMP];   
  Float_t vtxcomp_klmn_tk_chi2[MAX_VTXCOMP][4];  

  Bool_t  vtxcomp_adpt_isvalid[MAX_VTXCOMP];  
  Float_t vtxcomp_adpt_x2dof[MAX_VTXCOMP];    
  Int_t   vtxcomp_adpt_ndof[MAX_VTXCOMP];     
  Float_t vtxcomp_adpt_x2prob[MAX_VTXCOMP];   
  Float_t vtxcomp_adpt_tk_chi2[MAX_VTXCOMP][4];  
  Float_t vtxcomp_adpt_tk_weight[MAX_VTXCOMP][4];

 private:
  GlobeCuts *gCUT;
  int debug_level;
  edm::InputTag trackColl;
};

#endif
