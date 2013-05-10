#ifndef GLOBEPHOTONS_H
#define GLOBEPHOTONS_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"

#include "RecoEgamma/EgammaTools/interface/ggPFPhotons.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

typedef math::XYZTLorentzVector LorentzVector;

class GlobePhotons {
 public:

  GlobePhotons(const edm::ParameterSet&, const char* n = "unused");
  ~GlobePhotons();

  void checkSetup(const edm::EventSetup&);
  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  int PhotonID(reco::PhotonRef, int, reco::VertexRef, bool, int a = -1);

  Float_t SumTrackPtInConeHgg(reco::PhotonRef, reco::VertexRef, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t WorstSumTrackPtInConeHgg(reco::PhotonRef, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t DeltaRToTrackHgg(reco::PhotonRef, edm::Handle<reco::GsfElectronCollection>, int);
  LorentzVector VertexCorrectedP4Hgg(reco::PhotonRef, reco::VertexRef);

  int photonCutLevel6(float, float, float, float, float, float, float, float);
  int photonCutLevel4(float, float, float, float, float, float, float, float);

  void setPhotonIDThresholds(const edm::ParameterSet&);

  //std::vector<float> pfTkIsoWithVertex(reco::PhotonRef, const reco::PFCandidateCollection*, float, float, std::vector<reco::PFCandidate::ParticleType>); 
  //float pfEcalIso(reco::PhotonRef, const reco::PFCandidateCollection*, float, float, float, float, float, float, std::vector<reco::PFCandidate::ParticleType>);
  //float pfHcalIso(reco::PhotonRef, const reco::PFCandidateCollection*, float, float, std::vector<reco::PFCandidate::ParticleType>);
  //bool sameParticle(const reco::PFCandidate&, const reco::PFCandidate&) const;  

  std::map<DetId, EcalRecHit> rechits_map_;

    // variables

  Int_t pho_n;
  //correction schemes
  Float_t pho_feta[MAX_PHOTONS][5];
  Float_t pho_crackcorr[MAX_PHOTONS];
  Float_t pho_localcorr[MAX_PHOTONS];

  Int_t pho_isPFPhoton[MAX_PHOTONS];
  Int_t pho_isPFElectron[MAX_PHOTONS];

  //fiducial flags
  Int_t pho_isEB[MAX_PHOTONS];
  Int_t pho_isEE[MAX_PHOTONS];
  Int_t pho_isEBGap[MAX_PHOTONS];
  Int_t pho_isEEGap[MAX_PHOTONS];
  Int_t pho_isEBEEGap[MAX_PHOTONS];
  Int_t pho_isEBEtaGap[MAX_PHOTONS];
  Int_t pho_isEBPhiGap[MAX_PHOTONS];
  Int_t pho_isEEDeeGap[MAX_PHOTONS];
  Int_t pho_isEERingGap[MAX_PHOTONS];

  //shower shape variables
  Float_t pho_see[MAX_PHOTONS];
  Float_t pho_sieie[MAX_PHOTONS];
  Float_t pho_sipip[MAX_PHOTONS];
  Float_t pho_sieip[MAX_PHOTONS];
  Float_t pho_e1x5[MAX_PHOTONS];
  Float_t pho_e1x3[MAX_PHOTONS];
  Float_t pho_e2x2[MAX_PHOTONS];
  Float_t pho_e3x3[MAX_PHOTONS];
  Float_t pho_e5x5[MAX_PHOTONS];
  Float_t pho_emaxxtal[MAX_PHOTONS];
  Float_t pho_hoe[MAX_PHOTONS];
  Float_t pho_h1oe[MAX_PHOTONS];
  Float_t pho_h2oe[MAX_PHOTONS];
  Float_t pho_hoe_bc[MAX_PHOTONS];
  Float_t pho_h1oe_bc[MAX_PHOTONS];
  Float_t pho_h2oe_bc[MAX_PHOTONS];
  Float_t pho_r1x5[MAX_PHOTONS];
  Float_t pho_r2x5[MAX_PHOTONS];
  Float_t pho_r9[MAX_PHOTONS];
  Float_t pho_zernike20[MAX_PHOTONS];
  Float_t pho_zernike42[MAX_PHOTONS];
  Float_t pho_eseffsixix[MAX_PHOTONS];
  Float_t pho_eseffsiyiy[MAX_PHOTONS];
  Float_t pho_r19[MAX_PHOTONS];
  Float_t pho_maxoraw[MAX_PHOTONS];
  Float_t pho_cep[MAX_PHOTONS];
  Float_t pho_lambdaratio[MAX_PHOTONS];
  Float_t pho_lambdadivcov[MAX_PHOTONS];
  Float_t pho_cep_global[MAX_PHOTONS];
  Float_t pho_lambdaratio_global[MAX_PHOTONS];
  Float_t pho_lambdadivcov_global[MAX_PHOTONS];
  Float_t pho_etawidth[MAX_PHOTONS];
  Float_t pho_brem[MAX_PHOTONS];
  Float_t pho_smaj[MAX_PHOTONS];
  Float_t pho_hcalbcsumetconedr03[MAX_PHOTONS];
  Float_t pho_hcalbc1sumetconedr03[MAX_PHOTONS];
  Float_t pho_hcalbc2sumetconedr03[MAX_PHOTONS];
  Float_t pho_hcalbcsumetconedr04[MAX_PHOTONS];
  Float_t pho_hcalbc1sumetconedr04[MAX_PHOTONS];
  Float_t pho_hcalbc2sumetconedr04[MAX_PHOTONS];

  Float_t pho_sc_time[MAX_PHOTONS];

// pi0 disc variable
  Float_t pho_pi0disc[MAX_PHOTONS];

// OutIn Conv trks variable
  Int_t pho_IsConvOutIn[MAX_PHOTONS];

  Float_t pho_must[MAX_PHOTONS];
  Int_t pho_mustnc[MAX_PHOTONS];
  Float_t pho_e2nd[MAX_PHOTONS];
  Float_t pho_e2x5right[MAX_PHOTONS];
  Float_t pho_e2x5left[MAX_PHOTONS];
  Float_t pho_e2x5top[MAX_PHOTONS];
  Float_t pho_e2x5max[MAX_PHOTONS];
  Float_t pho_e2x5bottom[MAX_PHOTONS];
  Float_t pho_eright[MAX_PHOTONS];
  Float_t pho_eleft[MAX_PHOTONS];
  Float_t pho_etop[MAX_PHOTONS];
  Float_t pho_ebottom[MAX_PHOTONS];

  Int_t pho_bieta[MAX_PHOTONS];
  Int_t pho_biphi[MAX_PHOTONS];
  Float_t pho_betacry[MAX_PHOTONS];
  Float_t pho_bphicry[MAX_PHOTONS];
  Float_t pho_bthetatilt[MAX_PHOTONS];
  Float_t pho_bphitilt[MAX_PHOTONS];

  Float_t pho_pfconvVtxZ[MAX_PHOTONS];
  Float_t pho_pfconvVtxZErr[MAX_PHOTONS];
  Int_t pho_pfMatch[MAX_PHOTONS];
  Int_t pho_hasConvPf[MAX_PHOTONS];
  Int_t pho_hasSLConvPf[MAX_PHOTONS];
  Int_t pho_PfEleVeto[MAX_PHOTONS];    
  Float_t pho_pfpresh1[MAX_PHOTONS];
  Float_t pho_pfpresh2[MAX_PHOTONS];
  Float_t pho_mustenergy[MAX_PHOTONS];
  Float_t pho_mustenergyout[MAX_PHOTONS];
  Float_t pho_mustEtout[MAX_PHOTONS];
  Float_t pho_pflowE[MAX_PHOTONS];
  Float_t pho_pfdeta[MAX_PHOTONS];
  Float_t pho_pfdphi[MAX_PHOTONS];
  Float_t pho_pfclusrms[MAX_PHOTONS];
  Float_t pho_pfclusrmsmust[MAX_PHOTONS];
  Float_t pho_pfClusECorr[MAX_PHOTONS];

  Float_t pho_pfRawEnergy[MAX_PHOTONS];
  Float_t pho_pfe2x2[MAX_PHOTONS];
  Float_t pho_pfe3x3[MAX_PHOTONS];
  Float_t pho_pfe5x5[MAX_PHOTONS];
  Float_t pho_pfsieie[MAX_PHOTONS];
  Float_t pho_pfsieip[MAX_PHOTONS];
  Float_t pho_pfsipip[MAX_PHOTONS];
  Float_t pho_pfemaxxtal[MAX_PHOTONS];
  Float_t pho_pfe2nd[MAX_PHOTONS];


  Float_t pho_e2overe9[MAX_PHOTONS];
  Float_t pho_seed_time[MAX_PHOTONS];
  Float_t pho_seed_outoftimechi2[MAX_PHOTONS];
  Float_t pho_seed_chi2[MAX_PHOTONS];
  Float_t pho_seed_recoflag[MAX_PHOTONS];
  Float_t pho_seed_severity[MAX_PHOTONS];

  //isolation variables
  Float_t pho_pfiso_myneutral01[MAX_PHOTONS];
  Float_t pho_pfiso_myneutral02[MAX_PHOTONS];
  Float_t pho_pfiso_myneutral03[MAX_PHOTONS];
  Float_t pho_pfiso_myneutral04[MAX_PHOTONS];
  Float_t pho_pfiso_myneutral05[MAX_PHOTONS];
  Float_t pho_pfiso_myneutral06[MAX_PHOTONS];

  Float_t pho_pfiso_myphoton01[MAX_PHOTONS]; 
  Float_t pho_pfiso_myphoton02[MAX_PHOTONS]; 
  Float_t pho_pfiso_myphoton03[MAX_PHOTONS];  
  Float_t pho_pfiso_myphoton04[MAX_PHOTONS];
  Float_t pho_pfiso_myphoton05[MAX_PHOTONS]; 
  Float_t pho_pfiso_myphoton06[MAX_PHOTONS]; 

  /*
  Float_t pho_pfiso_egphoton01[MAX_PHOTONS]; 
  Float_t pho_pfiso_egphoton02[MAX_PHOTONS]; 
  Float_t pho_pfiso_egphoton03[MAX_PHOTONS];  
  Float_t pho_pfiso_egphoton04[MAX_PHOTONS];
  Float_t pho_pfiso_egphoton05[MAX_PHOTONS]; 
  Float_t pho_pfiso_egphoton06[MAX_PHOTONS]; 

  Float_t pho_pfiso_barephoton01[MAX_PHOTONS]; 
  Float_t pho_pfiso_barephoton02[MAX_PHOTONS]; 
  Float_t pho_pfiso_barephoton03[MAX_PHOTONS];  
  Float_t pho_pfiso_barephoton04[MAX_PHOTONS];
  Float_t pho_pfiso_barephoton05[MAX_PHOTONS]; 
  Float_t pho_pfiso_barephoton06[MAX_PHOTONS]; 
  */
  std::vector<std::vector<float> >* pho_pfiso_mycharged01;
  std::vector<std::vector<float> >* pho_pfiso_mycharged02;
  std::vector<std::vector<float> >* pho_pfiso_mycharged03;
  std::vector<std::vector<float> >* pho_pfiso_mycharged04;
  std::vector<std::vector<float> >* pho_pfiso_mycharged05;
  std::vector<std::vector<float> >* pho_pfiso_mycharged06;

  /*
  std::vector<std::vector<float> >* pho_pfiso_egcharged01;
  std::vector<std::vector<float> >* pho_pfiso_egcharged02;
  std::vector<std::vector<float> >* pho_pfiso_egcharged03;
  std::vector<std::vector<float> >* pho_pfiso_egcharged04;
  std::vector<std::vector<float> >* pho_pfiso_egcharged05;
  std::vector<std::vector<float> >* pho_pfiso_egcharged06;

  std::vector<std::vector<float> >* pho_pfiso_barecharged01;
  std::vector<std::vector<float> >* pho_pfiso_barecharged02;
  std::vector<std::vector<float> >* pho_pfiso_barecharged03;
  std::vector<std::vector<float> >* pho_pfiso_barecharged04;
  std::vector<std::vector<float> >* pho_pfiso_barecharged05;
  std::vector<std::vector<float> >* pho_pfiso_barecharged06;
  */
  //Float_t pho_pfiso_neutral03[MAX_PHOTONS];
  //Float_t pho_pfiso_charged03[MAX_PHOTONS];
  //Float_t pho_pfiso_photon03[MAX_PHOTONS];  
  //Float_t pho_pfiso_neutral04[MAX_PHOTONS];
  //Float_t pho_pfiso_charged04[MAX_PHOTONS];
  //Float_t pho_pfiso_photon04[MAX_PHOTONS];  

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
  //std::vector<std::vector<float> >* pho_frixiso;

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
  Float_t pho_regr_energy[MAX_PHOTONS];
  Float_t pho_regr_energyerr[MAX_PHOTONS];

  //Int_t pho_id_4cat[MAX_PHOTONS][100];  
  //Int_t pho_id_6cat[MAX_PHOTONS][100];
  //Int_t pho_id_6catpf[MAX_PHOTONS][100];
  
  std::vector<std::vector<UInt_t> >* pho_schits;
  std::vector<std::vector<UInt_t> >* pho_bchits;

  TClonesArray *pho_p4;
  TClonesArray *pho_calopos;

  TClonesArray *pho_conv_vtx;
  TClonesArray *pho_conv_pair_momentum;
  TClonesArray *pho_conv_refitted_momentum;
  TClonesArray *pho_conv_vertexcorrected_p4;

 private:
  GBRForest* PFLCBarrel;
  GBRForest* PFLCEndcap; 
  bool isInitialized;
  const char* nome;
  GlobeCuts *gCUT;
  GlobeEcalClusters *gES;
  edm::InputTag photonCollStd;
  edm::InputTag photonCollPf;
  CiCPhotonID* cicPhotonId;

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
  edm::InputTag ecalHitESColl;
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
  edm::InputTag pfColl;
  std::vector<edm::InputTag> inputTagIsoVals03_;
  std::vector<edm::InputTag> inputTagIsoVals04_;

  int debug_level;
  bool doFastSim;
  bool doAodSim;
  double rho; 
  
  edm::Handle<reco::VertexCollection> hVertex;
  edm::Handle<reco::TrackCollection> tkHandle;
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  edm::Handle<reco::PFCandidateCollection> pfCollection;
  edm::ESHandle<CaloGeometry> theCaloGeom_;

  const HBHERecHitCollection* hithbhe_;
// Correction Schemes
  EcalClusterFunctionBaseClass *fEtaCorr;
  EcalClusterFunctionBaseClass *CrackCorr;
  EcalClusterFunctionBaseClass *LocalCorr;

  /*
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
  */

  CaloGeometry geometry;
  const EcalSeverityLevelAlgo* sevLevel; 
  edm::ESHandle<TransientTrackBuilder> theTTkBuilder;

  ElectronHcalHelper *hcalHelper, *hcalHelperPflow;
  ElectronHcalHelper::Configuration hcalCfg, hcalCfgPflow;

  EGEnergyCorrector ecorr_;
  
  bool energyCorrectionsFromDB;
  std::string energyRegFilename;
  std::string regressionVersion;
};

#endif
