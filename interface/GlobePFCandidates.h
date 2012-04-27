#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEPFCANDIDATES_H
#define GLOBEPFCANDIDATES_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobePFCandidates {
 public:
  
  GlobePFCandidates(const edm::ParameterSet&);
  virtual ~GlobePFCandidates() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&, GlobeTracks*, GlobeMuons*, GlobePhotons*);

  // variables
  Int_t pfcand_n;
  Int_t pfcand_pdgid[MAX_PFCANDS];
  Int_t pfcand_tkind[MAX_PFCANDS];
  Int_t pfcand_gsfind[MAX_PFCANDS];
  Int_t pfcand_muind[MAX_PFCANDS];
  Float_t pfcand_ecalEnergy[MAX_PFCANDS];
  Float_t pfcand_hcalEnergy[MAX_PFCANDS];
  Float_t pfcand_rawEcalEnergy[MAX_PFCANDS];
  Float_t pfcand_rawHcalEnergy[MAX_PFCANDS];
  Float_t pfcand_ps1Energy[MAX_PFCANDS];
  Float_t pfcand_ps2Energy[MAX_PFCANDS];
  Float_t pfcand_momErr[MAX_PFCANDS];
  Float_t pfcand_mva_e_pi[MAX_PFCANDS];
  Float_t pfcand_mva_e_mu[MAX_PFCANDS];
  Float_t pfcand_mva_pi_mu[MAX_PFCANDS];
  Float_t pfcand_mva_nothing_gamma[MAX_PFCANDS];
  Float_t pfcand_mva_nothing_nh[MAX_PFCANDS];
  Float_t pfcand_mva_gamma_nh[MAX_PFCANDS];
  Float_t pfcand_vz[MAX_PFCANDS];

  Int_t pho_n;
  //Int_t pfcandtimespho_n;
  UInt_t pfcand_overlappho[MAX_PFCANDS];
  UInt_t pfcand_ispu[MAX_PFCANDS];
  //Int_t pfcand_overlappho[MAX_PFCANDS][MAX_PHOTONS];

  TClonesArray *pfcand_p4;
  TClonesArray *pfcand_poscalo;
  TClonesArray *pfcand_posvtx;

 private:
  GlobeCuts *gCUT;
  edm::InputTag pfColl;
  edm::InputTag photonCollStd;
  edm::InputTag electronCollStd;
  double PFIsoOuterConeSize;
  int debug_level;
};

#endif
