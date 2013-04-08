#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEMUONS_H
#define GLOBEMUONS_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeMuons {
 public:
  
  GlobeMuons(const edm::ParameterSet&, const char* n = "unused");
  virtual ~GlobeMuons();

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t mu_n;

  Float_t mu_em[MAX_MUONS];
  Float_t mu_had[MAX_MUONS];
  Float_t mu_ho[MAX_MUONS];
  Float_t mu_emS9[MAX_MUONS];
  Float_t mu_hadS9[MAX_MUONS];
  Float_t mu_hoS9[MAX_MUONS];
  Float_t mu_dz[MAX_MUONS];
  Float_t mu_d0[MAX_MUONS];
  Float_t mu_dxy[MAX_MUONS];
  Float_t mu_dzerr[MAX_MUONS];
  Float_t mu_d0err[MAX_MUONS];
  Float_t mu_chi2[MAX_MUONS];
  Float_t mu_dof[MAX_MUONS];
  Float_t mu_tkpterr[MAX_MUONS];
  Float_t mu_ecaliso03[MAX_MUONS];
  Float_t mu_hcaliso03[MAX_MUONS];
  Float_t mu_tkiso03[MAX_MUONS];
  Float_t mu_D0Vtx[MAX_MUONS][100];
  Float_t mu_DZVtx[MAX_MUONS][100];


  Int_t mu_charge[MAX_MUONS];
  Int_t mu_losthits[MAX_MUONS];
  Int_t mu_validhits[MAX_MUONS];
  Int_t mu_tkind[MAX_MUONS];
  Int_t mu_nmatches[MAX_MUONS];
  Int_t mu_staind[MAX_MUONS];
  Int_t mu_type[MAX_MUONS];
  Int_t mu_innerhits[MAX_MUONS];
  Int_t mu_pixelhits[MAX_MUONS];
  Int_t mu_validChmbhits[MAX_MUONS];

  Float_t mu_iso[MAX_MUONS];

  TClonesArray *mu_p4;
  TClonesArray *mu_momvtx;
  TClonesArray *mu_posvtx;
  TClonesArray *mu_poshcal;
  TClonesArray *mu_posecal;
  
 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag muonColl, trackColl, staMuonColl, vertexColl;
  edm::ParameterSet psetTAP;

  TrackAssociatorParameters* theAssociatorParameters;
  TrackDetectorAssociator* theAssociator;
  int debug_level;
  bool doAodSim;
  bool doEgammaSummer09Skim;
};

#endif
