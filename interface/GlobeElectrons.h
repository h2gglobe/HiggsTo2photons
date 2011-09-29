#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEELECTRONS_H
#define GLOBEELECTRONS_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeElectrons {
 public:
  
  GlobeElectrons(const edm::ParameterSet&, const char* n = "std");
  virtual ~GlobeElectrons() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  bool analyze_pf(const edm::Event&, const edm::EventSetup&);
  void initialize_branches(int electron_number);
  std::map<DetId, EcalRecHit> rechits_map_;

  std::pair<unsigned int, float> sharedHits(const reco::Track& trackA, const reco::Track& trackB);

  // variables
  Int_t el_n;
  Float_t el_enearbcopin[MAX_ELECTRONS];
  Float_t el_eseedopout[MAX_ELECTRONS];
  Float_t el_eopin[MAX_ELECTRONS];
  Float_t el_pout[MAX_ELECTRONS];
  Float_t el_pin[MAX_ELECTRONS];
  Float_t el_fbrem[MAX_ELECTRONS];
  Int_t el_nbrem[MAX_ELECTRONS];
  Int_t el_1pxb[MAX_ELECTRONS];
  Int_t el_1pxf[MAX_ELECTRONS];

  Float_t el_hoe[MAX_ELECTRONS];
  Float_t el_hoed1[MAX_ELECTRONS];
  Float_t el_hoed2[MAX_ELECTRONS];
  Float_t el_h[MAX_ELECTRONS];

  Float_t el_detain[MAX_ELECTRONS];
  Float_t el_dphiin[MAX_ELECTRONS];
  Float_t el_detaout[MAX_ELECTRONS];
  Float_t el_dphiout[MAX_ELECTRONS];
  Float_t el_z0[MAX_ELECTRONS];
  Float_t el_d0[MAX_ELECTRONS];
  Float_t el_z0err[MAX_ELECTRONS];
  Float_t el_d0err[MAX_ELECTRONS];
  Float_t el_chi2[MAX_ELECTRONS];
  Float_t el_dof[MAX_ELECTRONS];
  Float_t el_e1x5[MAX_ELECTRONS];
  Float_t el_e5x5[MAX_ELECTRONS];
  Float_t el_sipip[MAX_ELECTRONS];
  Float_t el_sieie[MAX_ELECTRONS];
  Float_t el_sieiesc[MAX_ELECTRONS];
  Float_t el_eseffsixix[MAX_ELECTRONS];
  Float_t el_eseffsiyiy[MAX_ELECTRONS];
  Float_t el_e2x5[MAX_ELECTRONS];
  Float_t el_esc[MAX_ELECTRONS];
  Float_t el_eseed[MAX_ELECTRONS];
  Float_t el_eseedopin[MAX_ELECTRONS];
  Float_t el_qoverperr[MAX_ELECTRONS];
  Float_t el_pterr[MAX_ELECTRONS];
  Float_t el_etaerr[MAX_ELECTRONS];
  Float_t el_phierr[MAX_ELECTRONS];
  Int_t el_class[MAX_ELECTRONS];
  Int_t el_charge[MAX_ELECTRONS];
  Int_t el_ch_gsf[MAX_ELECTRONS];
  Int_t el_ch_scpix[MAX_ELECTRONS];
  Int_t el_losthits[MAX_ELECTRONS];
  Int_t el_validhits[MAX_ELECTRONS];
  Int_t el_hp_expin[MAX_ELECTRONS];
  Int_t el_hp_expout[MAX_ELECTRONS];

  Int_t el_scind[MAX_ELECTRONS];
  Int_t el_crack[MAX_ELECTRONS];
  Int_t el_tkind[MAX_ELECTRONS];
  Int_t el_conv[MAX_ELECTRONS];
  
  Int_t el_nambtk[MAX_ELECTRONS];

  Float_t el_pfiso_myneutral03[MAX_ELECTRONS];
  Float_t el_pfiso_mycharged03[MAX_ELECTRONS];
  Float_t el_pfiso_myphoton03[MAX_ELECTRONS];
  Float_t el_pfiso_myneutral04[MAX_ELECTRONS];
  Float_t el_pfiso_mycharged04[MAX_ELECTRONS];
  Float_t el_pfiso_myphoton04[MAX_ELECTRONS];
  Float_t el_pfiso_neutral[MAX_ELECTRONS];
  Float_t el_pfiso_charged[MAX_ELECTRONS];
  Float_t el_pfiso_photon[MAX_ELECTRONS];
  Float_t el_tkiso03[MAX_ELECTRONS];
  Float_t el_ecaliso03[MAX_ELECTRONS];
  Float_t el_hcaliso03[MAX_ELECTRONS];
  Float_t el_tkiso04[MAX_ELECTRONS];
  Float_t el_ecaliso04[MAX_ELECTRONS];
  Float_t el_hcaliso04[MAX_ELECTRONS];

  Float_t el_mva[MAX_ELECTRONS];
  Bool_t el_ecaldrv[MAX_ELECTRONS];
  Bool_t el_tkdrv[MAX_ELECTRONS];
  Float_t el_ip_ctf[MAX_ELECTRONS];
  Float_t el_ip_gsf[MAX_ELECTRONS];
  Float_t el_dist[MAX_ELECTRONS];
  Float_t el_dcot[MAX_ELECTRONS];

  std::vector<std::vector<int> >* el_catbased;

  TClonesArray *el_sc;
  TClonesArray *el_p4;
  TClonesArray *el_momvtx;
  TClonesArray *el_momvtxconst;
  TClonesArray *el_momcalo;
  TClonesArray *el_momout;
  TClonesArray *el_posvtx;
  TClonesArray *el_poscalo;

 private:
  const char* nome; 
  bool doFastSim;
  bool doAodSim;
  GlobeCuts *gCUT;
  edm::InputTag electronColl, trackColl, trackColl2, vertexColl, beamSpotColl, conversionColl, pfColl;
  GlobeEcalClusters *gES;
  std::vector<edm::InputTag> eIDLabels;
  int debug_level;

  // SUPER CLUSTERS
  edm::InputTag hybridSuperClusterColl;
  edm::InputTag endcapSuperClusterColl;  
  edm::InputTag ecalHitEBColl;
  edm::InputTag ecalHitEEColl;
  edm::InputTag ecalHitESColl;
  edm::InputTag hcalHitColl;
};

#endif
