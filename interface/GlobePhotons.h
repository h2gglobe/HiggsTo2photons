#ifndef GLOBEPHOTONS_H
#define GLOBEPHOTONS_H

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

class GlobePhotons {
 public:
  
  GlobePhotons(const edm::ParameterSet&, const char* n = "unused");
  virtual ~GlobePhotons() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  double getHoE(GlobalPoint, double, const edm::Event&, const edm::EventSetup&);

  // variables

  Int_t pho_n;
  //fiducial flags
  Int_t pho_isEB[MAX_PHOTONS];
  Int_t pho_isEE[MAX_PHOTONS];
  Int_t pho_isEBGap[MAX_PHOTONS];
  Int_t pho_isEEGap[MAX_PHOTONS];
  Int_t pho_isEBEEGap[MAX_PHOTONS];

  //shower shape variables
  Float_t pho_see[MAX_PHOTONS];
  Float_t pho_sieie[MAX_PHOTONS];
  Float_t pho_e1x5[MAX_PHOTONS];
  Float_t pho_e2x5[MAX_PHOTONS];
  Float_t pho_e3x3[MAX_PHOTONS];
  Float_t pho_e5x5[MAX_PHOTONS];
  Float_t pho_emaxxtal[MAX_PHOTONS];
  Float_t pho_hoe[MAX_PHOTONS];
  Float_t pho_h1oe[MAX_PHOTONS];
  Float_t pho_h2oe[MAX_PHOTONS];
  Float_t pho_r1x5[MAX_PHOTONS];
  Float_t pho_r2x5[MAX_PHOTONS];
  Float_t pho_r9[MAX_PHOTONS];

  //isolation variables
  Float_t pho_ecalsumetconedr04[MAX_PHOTONS];
  Float_t pho_hcalsumetconedr04[MAX_PHOTONS];
  Float_t pho_hcal1sumetconedr04[MAX_PHOTONS];
  Float_t pho_hcal2sumetconedr04[MAX_PHOTONS];
  Float_t pho_trksumptsolidconedr04[MAX_PHOTONS];
  Float_t pho_trksumpthollowconedr04[MAX_PHOTONS];
  Float_t pho_ntrksolidconedr04[MAX_PHOTONS];
  Float_t pho_ntrkhollowconedr04[MAX_PHOTONS];
  Float_t pho_ecalsumetconedr03[MAX_PHOTONS];
  Float_t pho_hcalsumetconedr03[MAX_PHOTONS];
  Float_t pho_hcal1sumetconedr03[MAX_PHOTONS];
  Float_t pho_hcal2sumetconedr03[MAX_PHOTONS];
  Float_t pho_trksumptsolidconedr03[MAX_PHOTONS];
  Float_t pho_trksumpthollowconedr03[MAX_PHOTONS];
  Float_t pho_ntrksolidconedr03[MAX_PHOTONS];
  Float_t pho_ntrkhollowconedr03[MAX_PHOTONS];

  Int_t pho_barrel[MAX_PHOTONS];
  Int_t pho_scind[MAX_PHOTONS];
  Int_t pho_haspixseed[MAX_PHOTONS];

  Int_t pho_hasconvtks[MAX_PHOTONS];
  Int_t pho_nconv[MAX_PHOTONS];
  Int_t pho_conv_ntracks[MAX_PHOTONS];
  Float_t pho_conv_pairinvmass[MAX_PHOTONS];
  Float_t pho_conv_paircotthetasep[MAX_PHOTONS];
  Float_t pho_conv_eoverp[MAX_PHOTONS];
  Float_t pho_conv_zofprimvtxfromtrks[MAX_PHOTONS];
  Float_t pho_conv_distofminapproach[MAX_PHOTONS];
  Float_t pho_conv_dphitrksatvtx[MAX_PHOTONS];
  Float_t pho_conv_dphitrksatecal[MAX_PHOTONS];
  Float_t pho_conv_detatrksatecal[MAX_PHOTONS];
  Float_t pho_conv_tk1_d0[MAX_PHOTONS];
  Float_t pho_conv_tk1_pout[MAX_PHOTONS];
  Float_t pho_conv_tk1_pin[MAX_PHOTONS];
  Float_t pho_conv_tk2_d0[MAX_PHOTONS];
  Float_t pho_conv_tk2_pout[MAX_PHOTONS];
  Float_t pho_conv_tk2_pin[MAX_PHOTONS];

  Float_t pho_conv_tk1_dz[MAX_PHOTONS];
  Float_t pho_conv_tk1_dzerr[MAX_PHOTONS];
  Int_t pho_conv_tk1_nh[MAX_PHOTONS];
  Float_t pho_conv_tk2_dz[MAX_PHOTONS];
  Float_t pho_conv_tk2_dzerr[MAX_PHOTONS];
  Int_t pho_conv_tk2_nh[MAX_PHOTONS];
  Int_t pho_conv_ch1ch2[MAX_PHOTONS];
  Float_t pho_conv_chi2[MAX_PHOTONS];
  Float_t pho_conv_chi2_probability[MAX_PHOTONS];
  Int_t pho_conv_validvtx[MAX_PHOTONS];
  Float_t pho_conv_MVALikelihood[MAX_PHOTONS];
  
  TClonesArray *pho_p4;
  TClonesArray *pho_calopos;

  TClonesArray *pho_conv_vtx;
  TClonesArray *pho_conv_pair_momentum;
  TClonesArray *pho_conv_refitted_momentum;
  TClonesArray *pho_conv_vertexcorrected_p4;

 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag photonCollStd;

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

  int debug_level;
  bool doConvertedPhoton;
  bool doFastSim;
  bool doAodSim;

  edm::ESHandle<CaloGeometry> theCaloGeom_;
  const HBHERecHitCollection* hithbhe_;
};

#endif
