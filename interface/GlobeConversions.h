#ifndef GLOBECONVERSIONS_H
#define GLOBECONVERSIONS_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollection.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Math/VectorUtil.h"
#include <iostream>

class GlobeAnalyzer;

class GlobeConversions {
 public:
  
  GlobeConversions(const edm::ParameterSet&, const char* n = "unused");
  virtual ~GlobeConversions() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables

  Int_t conv_n;
  TClonesArray *conv_p4;
  TClonesArray *conv_vtx;
  TClonesArray *conv_pair_momentum;
  TClonesArray *conv_refitted_momentum;
  TClonesArray *conv_singleleg_momentum;


  Int_t conv_barrel[MAX_CONVERTEDPHOTONS];
  Int_t conv_scind[MAX_CONVERTEDPHOTONS];
  Int_t conv_ntracks[MAX_CONVERTEDPHOTONS];
  Float_t conv_pairinvmass[MAX_CONVERTEDPHOTONS];
  Float_t conv_paircotthetasep[MAX_CONVERTEDPHOTONS];
  Float_t conv_eoverp[MAX_CONVERTEDPHOTONS];
  Float_t conv_distofminapproach[MAX_CONVERTEDPHOTONS];
  Float_t conv_dphitrksatvtx[MAX_CONVERTEDPHOTONS];
  Float_t conv_dphitrksatecal[MAX_CONVERTEDPHOTONS];
  Float_t conv_detatrksatecal[MAX_CONVERTEDPHOTONS];
  std::vector<std::vector<int> >* conv_quality;
  Int_t conv_type[MAX_CONVERTEDPHOTONS];
  Float_t conv_dxy[MAX_CONVERTEDPHOTONS];
  Float_t conv_dz[MAX_CONVERTEDPHOTONS];
  Float_t conv_lxy[MAX_CONVERTEDPHOTONS];
  Float_t conv_lz[MAX_CONVERTEDPHOTONS];
  std::vector<std::vector<unsigned short> >* conv_nHitsBeforeVtx;
  uint8_t conv_nSharedHits[MAX_CONVERTEDPHOTONS];
  Float_t conv_zofprimvtxfromtrks[MAX_CONVERTEDPHOTONS];

  /// track quantities
  Float_t conv_tk1_d0[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_pout[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_pin[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_d0[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_pout[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_pin[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_dz[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_dzerr[MAX_CONVERTEDPHOTONS];
  Short_t conv_tk1_nh[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_dz[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_dzerr[MAX_CONVERTEDPHOTONS];
  Short_t conv_tk2_nh[MAX_CONVERTEDPHOTONS];
  Int_t conv_ch1ch2[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_pterr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_pterr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_etaerr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_etaerr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_thetaerr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_thetaerr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_phierr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_phierr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk1_lambdaerr[MAX_CONVERTEDPHOTONS];
  Float_t conv_tk2_lambdaerr[MAX_CONVERTEDPHOTONS];

  /// vertex quanities
  Float_t conv_vtx_xErr[MAX_CONVERTEDPHOTONS];
  Float_t conv_vtx_yErr[MAX_CONVERTEDPHOTONS];
  Float_t conv_vtx_zErr[MAX_CONVERTEDPHOTONS];
  Float_t conv_chi2[MAX_CONVERTEDPHOTONS];
  Float_t conv_chi2_probability[MAX_CONVERTEDPHOTONS];
  Int_t conv_validvtx[MAX_CONVERTEDPHOTONS];
  Float_t conv_MVALikelihood[MAX_CONVERTEDPHOTONS];

  Float_t conv_vtxProb[MAX_CONVERTEDPHOTONS];
  Int_t conv_nHitsMax[MAX_CONVERTEDPHOTONS];
  Int_t conv_eleind[MAX_CONVERTEDPHOTONS];
  Float_t conv_lxy_bs[MAX_CONVERTEDPHOTONS];

 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag photonCollStd;
  edm::InputTag allConversionsColl;

  // SUPER CLUSTERS
  edm::InputTag hybridSuperClusterColl;
  edm::InputTag barrelSuperClusterColl;
  edm::InputTag endcapSuperClusterColl;
  // BASIC CLUSTERS
  edm::InputTag barrelHybridClusterShapeColl; 
  edm::InputTag endcapBasicClusterShapeColl; 

// ECAL HITS
  edm::InputTag ecalHitEBColl;
  edm::InputTag ecalHitEEColl;
// HCAL HITS
  edm::InputTag hcalBEColl;
  edm::InputTag hcalFColl;
  edm::InputTag hcalHoColl; 
// Particle Flow
  //edm::InputTag pfColl;
  edm::InputTag pfPhotonsColl;

  edm::InputTag beamSpotColl;
  edm::InputTag eleColl;

  int debug_level;
  bool doAodSim;

  edm::ESHandle<CaloGeometry> theCaloGeom_;
  const HBHERecHitCollection* hithbhe_;
};

#endif
