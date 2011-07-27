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

#include "DataFormats/VertexReco/interface/VertexFwd.h"

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

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterEnergyCorrection.h"
#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterCrackCorrection.h"
#include "RecoEcal/EgammaCoreTools/plugins/EcalClusterLocalContCorrection.h"


#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <iostream>

typedef math::XYZTLorentzVector LorentzVector;

class GlobePhotons {
 public:

  GlobePhotons(const edm::ParameterSet&, const char* n = "unused");
  virtual ~GlobePhotons() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  float hoeCalculator(const reco::BasicCluster*, const CaloGeometry&,
                      const edm::Event&, const edm::EventSetup&);

  int PhotonID(reco::PhotonRef, int, reco::VertexRef, bool, int a = -1);

  Float_t SumTrackPtInConeHgg(reco::PhotonRef, reco::VertexRef, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t WorstSumTrackPtInConeHgg(reco::PhotonRef, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t DeltaRToTrackHgg(reco::PhotonRef, edm::Handle<reco::GsfElectronCollection>, int);
  LorentzVector VertexCorrectedP4Hgg(reco::PhotonRef, reco::VertexRef);

  int photonCutLevel6(float, float, float, float, float, float, float, float);
  int photonCutLevel4(float, float, float, float, float, float, float, float);

  void setPhotonIDThresholds(const edm::ParameterSet&);

  // variables

  Int_t pho_n;
  //correction schemes
  Float_t pho_feta[MAX_PHOTONS][5];
  Float_t pho_crackcorr[MAX_PHOTONS];
  Float_t pho_localcorr[MAX_PHOTONS];

  //fiducial flags
  Int_t pho_isEB[MAX_PHOTONS];
  Int_t pho_isEE[MAX_PHOTONS];
  Int_t pho_isEBGap[MAX_PHOTONS];
  Int_t pho_isEEGap[MAX_PHOTONS];
  Int_t pho_isEBEEGap[MAX_PHOTONS];

  //shower shape variables
  Float_t pho_see[MAX_PHOTONS];
  Float_t pho_sieie[MAX_PHOTONS];
  Float_t pho_sipip[MAX_PHOTONS];
  Float_t pho_sieip[MAX_PHOTONS];
  Float_t pho_e1x5[MAX_PHOTONS];
  Float_t pho_e2x5[MAX_PHOTONS];
  Float_t pho_e3x3[MAX_PHOTONS];
  Float_t pho_e5x5[MAX_PHOTONS];
  Float_t pho_emaxxtal[MAX_PHOTONS];
  Float_t pho_hoe[MAX_PHOTONS];
  Float_t pho_h[MAX_PHOTONS];
  Float_t pho_h1oe[MAX_PHOTONS];
  Float_t pho_h2oe[MAX_PHOTONS];
  Float_t pho_r1x5[MAX_PHOTONS];
  Float_t pho_r2x5[MAX_PHOTONS];
  Float_t pho_r9[MAX_PHOTONS];
  Float_t pho_zernike20[MAX_PHOTONS];
  Float_t pho_zernike42[MAX_PHOTONS];

// pi0 disc variable
  Float_t pho_pi0disc[MAX_PHOTONS];

// OutIn Conv trks variable
  Int_t pho_IsConvOutIn[MAX_PHOTONS];

  Float_t pho_e2nd[MAX_PHOTONS];
  Float_t pho_e2x5right[MAX_PHOTONS];
  Float_t pho_e2x5left[MAX_PHOTONS];
  Float_t pho_e2x5Top[MAX_PHOTONS];
  Float_t pho_e2x5bottom[MAX_PHOTONS];
  Float_t pho_eright[MAX_PHOTONS];
  Float_t pho_eleft[MAX_PHOTONS];
  Float_t pho_etop[MAX_PHOTONS];
  Float_t pho_ebottom[MAX_PHOTONS];

  Float_t pho_e2overe9[MAX_PHOTONS];
  Float_t pho_seed_time[MAX_PHOTONS];
  Float_t pho_seed_outoftimechi2[MAX_PHOTONS];
  Float_t pho_seed_chi2[MAX_PHOTONS];
  Float_t pho_seed_recoflag[MAX_PHOTONS];
  Float_t pho_seed_severity[MAX_PHOTONS];

  //isolation variables
  Float_t pho_pfiso_neutral03[MAX_PHOTONS];
  Float_t pho_pfiso_charged03[MAX_PHOTONS];
  Float_t pho_pfiso_photon03[MAX_PHOTONS];
  Float_t pho_pfiso_neutral04[MAX_PHOTONS];
  Float_t pho_pfiso_charged04[MAX_PHOTONS];
  Float_t pho_pfiso_photon04[MAX_PHOTONS];
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
  Int_t pho_isconv[MAX_PHOTONS];

  Float_t pho_residCorrEnergy[MAX_PHOTONS];
  Float_t pho_residCorrResn[MAX_PHOTONS];

  Int_t pho_id[MAX_PHOTONS];

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

  edm::InputTag convertedPhotonColl;
  edm::InputTag beamSpotColl;
  edm::InputTag electronColl;

  edm::InputTag rhoCollection;
  edm::InputTag vtxCollection;
  edm::InputTag tkCollection;
  edm::InputTag hcalHitColl;
  std::vector<edm::InputTag> inputTagIsoVals03_;
  std::vector<edm::InputTag> inputTagIsoVals04_;

  int debug_level;
  bool doFastSim;
  bool doAodSim;
  double rho; 
  
  edm::Handle<reco::VertexCollection> hVertex;
  edm::Handle<reco::TrackCollection> tkHandle;
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  edm::ESHandle<CaloGeometry> theCaloGeom_;

  const HBHERecHitCollection* hithbhe_;
// Correction Schemes
  EcalClusterFunctionBaseClass *fEtaCorr;
  EcalClusterFunctionBaseClass *CrackCorr;
  EcalClusterFunctionBaseClass *LocalCorr;

  // PhotonID thresholds
  std::vector<double> cutsubleadisosumoet6c[12];
  std::vector<double> cutsubleadisosumoetbad6c[12];
  std::vector<double> cutsubleadtrkisooetom6c[12];
  std::vector<double> cutsubleadsieie6c[12];
  std::vector<double> cutsubleadhovere6c[12];
  std::vector<double> cutsubleadr96c[12];
  std::vector<double> cutsublead_drtotk6c[12];
  
  std::vector<double> cutsubleadisosumoet[12];
  std::vector<double> cutsubleadisosumoetbad[12];
  std::vector<double> cutsubleadtrkisooetom[12];
  std::vector<double> cutsubleadsieie[12];
  std::vector<double> cutsubleadhovere[12];
  std::vector<double> cutsubleadr9[12];
  std::vector<double> cutsublead_drtotk[12];
};

#endif
