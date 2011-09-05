#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

///////////////
// pi0 disc
#include "DataFormats/EgammaCandidates/interface/PhotonPi0DiscriminatorAssociation.h"
// OutIn Conv tracks
//
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionTrackEcalImpactPoint.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionTrackPairFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionVertexFinder.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConvertedPhotonProducer.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
//
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/EgammaTrackReco/interface/TrackCaloClusterAssociation.h"
/////////////////

#include "HiggsAnalysis/HiggsTo2photons/interface/PhotonFixCMS.h"
#include "DataFormats/Math/interface/deltaR.h"


GlobePhotons::GlobePhotons(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  debug_level = iConfig.getParameter<int>("Debug_Level");
  doFastSim = iConfig.getParameter<bool>("doFastSim");
  doAodSim = iConfig.getParameter<bool>("doAodSim");

  // PHOTONS 
  photonCollStd =  iConfig.getParameter<edm::InputTag>("PhotonCollStd");

  // SUPER CLUSTERS
  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");

  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");

  hcalBEColl =  iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");
  hcalFColl =  iConfig.getParameter<edm::InputTag>("HcalHitsFColl");
  hcalHoColl =  iConfig.getParameter<edm::InputTag>("HcalHitsHoColl");

  convertedPhotonColl =  iConfig.getParameter<edm::InputTag>("ConvertedPhotonColl");
  beamSpotColl =  iConfig.getParameter<edm::InputTag>("BeamSpot");
  electronColl =  iConfig.getParameter<edm::InputTag>("ElectronColl_std");

  rhoCollection = iConfig.getParameter<edm::InputTag>("rhoCorrection");
  vtxCollection = iConfig.getParameter<edm::InputTag>("VertexColl_std");
  tkCollection  = iConfig.getParameter<edm::InputTag>("tkColl");

  hcalHitColl = iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");
  pfColl = iConfig.getParameter<edm::InputTag>("PFCandidateColl");

  edm::ParameterSet isoVals03  = iConfig.getParameter<edm::ParameterSet> ("isolationValues03");
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfChargedHadrons"));
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfPhotons"));
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfNeutralHadrons"));
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfPhotonsNoveto"));

  edm::ParameterSet isoVals04  = iConfig.getParameter<edm::ParameterSet> ("isolationValues04");
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfChargedHadrons"));
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfPhotons"));
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfNeutralHadrons"));
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfPhotonsNoveto"));

  // get the Correction Functions
  fEtaCorr  = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrection",iConfig);
  CrackCorr = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection",iConfig);
  LocalCorr = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);

  // set Hgg PhotonID thresholds
  setPhotonIDThresholds(iConfig);

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);

  pho_pfiso_mycharged03 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged04 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged03_noveto = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged04_noveto = new std::vector<std::vector<float> >();

}


void GlobePhotons::setPhotonIDThresholds(const edm::ParameterSet& iConfig) {

  const edm::ParameterSet gammaIDCuts = iConfig.getParameter<edm::ParameterSet>("hggPhotonIDConfiguration") ;
  char a[100];

  for (int lev = 0; lev < 11; ++lev) {
    sprintf(a, "cutsubleadisosumoet6c%d", lev);
    cutsubleadisosumoet6c[lev]     = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadisosumoetbad6c%d", lev);
    cutsubleadisosumoetbad6c[lev]  = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadtrkisooetom6c%d", lev);
    cutsubleadtrkisooetom6c[lev]   = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadsieie6c%d", lev);
    cutsubleadsieie6c[lev]         = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadhovere6c%d", lev);
    cutsubleadhovere6c[lev]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadr96c%d", lev);
    cutsubleadr96c[lev]            = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsublead_drtotk_25_996c%d", lev);
    cutsublead_drtotk6c[lev] = gammaIDCuts.getParameter<std::vector<double> >(a);
    
    sprintf(a, "cutsubleadisosumoet%d", lev);
    cutsubleadisosumoet[lev]     = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadisosumoetbad%d", lev);
    cutsubleadisosumoetbad[lev]  = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadtrkisooetom%d", lev);
    cutsubleadtrkisooetom[lev]   = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadsieie%d", lev);
    cutsubleadsieie[lev]         = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadhovere%d", lev);
    cutsubleadhovere[lev]        = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsubleadr9%d", lev);
    cutsubleadr9[lev]            = gammaIDCuts.getParameter<std::vector<double> >(a);
    sprintf(a, "cutsublead_drtotk_25_99%d", lev);
    cutsublead_drtotk[lev] = gammaIDCuts.getParameter<std::vector<double> >(a);
  }
}


void GlobePhotons::defineBranch(TTree* tree) {

  pho_p4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
  pho_calopos = new TClonesArray("TVector3", MAX_PHOTONS);

  tree->Branch("pho_n", &pho_n, "pho_n/I");

  //Correction Schemes
  tree->Branch("pho_feta",&pho_feta,"pho_feta[pho_n][5]/F");
  tree->Branch("pho_crackcorr",&pho_crackcorr,"pho_crackcorr[pho_n]/F");
  tree->Branch("pho_localcorr",&pho_localcorr,"pho_localcorr[pho_n]/F");

  //fiducial flags
  tree->Branch("pho_isEB",&pho_isEB,"pho_isEB[pho_n]/I");
  tree->Branch("pho_isEE",&pho_isEE,"pho_isEE[pho_n]/I");
  tree->Branch("pho_isEBGap",&pho_isEBGap,"pho_isEBGap[pho_n]/I");
  tree->Branch("pho_isEEGap",&pho_isEEGap,"pho_isEEGap[pho_n]/I");
  tree->Branch("pho_isEBEEGap",&pho_isEBEEGap,"pho_isEBEEGap[pho_n]/I");

  //shower shape variables
  tree->Branch("pho_see",&pho_see,"pho_see[pho_n]/F");
  tree->Branch("pho_sieie",&pho_sieie,"pho_sieie[pho_n]/F");
  tree->Branch("pho_e1x5",&pho_e1x5,"pho_e1x5[pho_n]/F");
  tree->Branch("pho_e2x5",&pho_e2x5,"pho_e2x5[pho_n]/F");
  tree->Branch("pho_e3x3",&pho_e3x3,"pho_e3x3[pho_n]/F");
  tree->Branch("pho_e5x5",&pho_e5x5,"pho_e5x5[pho_n]/F");
  tree->Branch("pho_emaxxtal",&pho_emaxxtal,"pho_emaxxtal[pho_n]/F");
  tree->Branch("pho_hoe",&pho_hoe,"pho_hoe[pho_n]/F");
  tree->Branch("pho_h", &pho_h,"pho_h[pho_n]/F");
  tree->Branch("pho_h1oe",&pho_h1oe,"pho_h1oe[pho_n]/F");
  tree->Branch("pho_h2oe",&pho_h2oe,"pho_h2oe[pho_n]/F");
  tree->Branch("pho_r1x5",&pho_r1x5,"pho_r1x5[pho_n]/F");
  tree->Branch("pho_r2x5",&pho_r2x5,"pho_r2x5[pho_n]/F");
  tree->Branch("pho_r9",&pho_r9,"pho_r9[pho_n]/F");
  
  // added by Aris
  // pi0 disc
  tree->Branch("pho_pi0disc",&pho_pi0disc,"pho_pi0disc[pho_n]/F");
  // OutIn Conv trks
  tree->Branch("pho_IsConvOutIn",&pho_IsConvOutIn,"pho_IsConvOutIn[pho_n]/I");
  ////////

  //isolation variables
  //tree->Branch("pho_pfiso_charged03", &pho_pfiso_charged03, "pho_pfiso_charged03[pho_n]/F");
  tree->Branch("pho_pfiso_neutral03", &pho_pfiso_neutral03, "pho_pfiso_neutral03[pho_n]/F");
  tree->Branch("pho_pfiso_photon03", &pho_pfiso_photon03, "pho_pfiso_photon03[pho_n]/F");  
  tree->Branch("pho_pfiso_photon03_noveto", &pho_pfiso_photon03_noveto, "pho_pfiso_photon03_noveto[pho_n]/F");

  //tree->Branch("pho_pfiso_charged04", &pho_pfiso_charged04, "pho_pfiso_charged04[pho_n]/F");
  tree->Branch("pho_pfiso_neutral04", &pho_pfiso_neutral04, "pho_pfiso_neutral04[pho_n]/F");
  tree->Branch("pho_pfiso_photon04", &pho_pfiso_photon04, "pho_pfiso_photon04[pho_n]/F");
  tree->Branch("pho_pfiso_photon04_noveto", &pho_pfiso_photon04_noveto, "pho_pfiso_photon04_noveto[pho_n]/F");

  tree->Branch("pho_pfiso_mycharged03", "std::vector<std::vector<float> >", &pho_pfiso_mycharged03);
  tree->Branch("pho_pfiso_mycharged04", "std::vector<std::vector<float> >", &pho_pfiso_mycharged04);
  tree->Branch("pho_pfiso_mycharged03_noveto", "std::vector<std::vector<float> >", &pho_pfiso_mycharged03_noveto);
  tree->Branch("pho_pfiso_mycharged04_noveto", "std::vector<std::vector<float> >", &pho_pfiso_mycharged04_noveto);

  tree->Branch("pho_ecalsumetconedr04",&pho_ecalsumetconedr04,"pho_ecalsumetconedr04[pho_n]/F");
  tree->Branch("pho_hcalsumetconedr04",&pho_hcalsumetconedr04,"pho_hcalsumetconedr04[pho_n]/F");
  tree->Branch("pho_hcal1sumetconedr04",&pho_hcal1sumetconedr04,"pho_hcal1sumetconedr04[pho_n]/F");
  tree->Branch("pho_hcal2sumetconedr04",&pho_hcal2sumetconedr04,"pho_hcal2sumetconedr04[pho_n]/F");
  tree->Branch("pho_trksumptsolidconedr04",&pho_trksumptsolidconedr04,"pho_trksumptsolidconedr04[pho_n]/F");
  tree->Branch("pho_trksumpthollowconedr04",&pho_trksumpthollowconedr04,"pho_trksumpthollowconedr04[pho_n]/F");
  tree->Branch("pho_ntrksolidconedr04",&pho_ntrksolidconedr04,"pho_ntrksolidconedr04[pho_n]/F");
  tree->Branch("pho_ntrkhollowconedr04",&pho_ntrkhollowconedr04,"pho_ntrkhollowconedr04[pho_n]/F");
  tree->Branch("pho_ecalsumetconedr03",&pho_ecalsumetconedr03,"pho_ecalsumetconedr03[pho_n]/F");
  tree->Branch("pho_hcalsumetconedr03",&pho_hcalsumetconedr03,"pho_hcalsumetconedr03[pho_n]/F");
  tree->Branch("pho_hcal1sumetconedr03",&pho_hcal1sumetconedr03,"pho_hcal1sumetconedr03[pho_n]/F");
  tree->Branch("pho_hcal2sumetconedr03",&pho_hcal2sumetconedr03,"pho_hcal2sumetconedr03[pho_n]/F");
  tree->Branch("pho_trksumptsolidconedr03",&pho_trksumptsolidconedr03,"pho_trksumptsolidconedr03[pho_n]/F");
  tree->Branch("pho_trksumpthollowconedr03",&pho_trksumpthollowconedr03,"pho_trksumpthollowconedr03[pho_n]/F");
  tree->Branch("pho_ntrksolidconedr03",&pho_ntrksolidconedr03,"pho_ntrksolidconedr03[pho_n]/F");
  tree->Branch("pho_ntrkhollowconedr03",&pho_ntrkhollowconedr03,"pho_ntrkhollowconedr03[pho_n]/F");

  tree->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);
  tree->Branch("pho_calopos", "TClonesArray", &pho_calopos, 32000, 0);
  tree->Branch("pho_barrel", &pho_barrel, "pho_barrel[pho_n]/I");
  tree->Branch("pho_scind", &pho_scind, "pho_scind[pho_n]/I");
  
  tree->Branch("pho_haspixseed",&pho_haspixseed,"pho_haspixseed[pho_n]/I");
  
  tree->Branch("pho_hasconvtks",&pho_hasconvtks,"pho_hasconvtks[pho_n]/I");
  tree->Branch("pho_nconv",&pho_nconv,"pho_nconv[pho_n]/I");
  tree->Branch("pho_conv_ntracks",&pho_conv_ntracks,"pho_conv_ntracks[pho_n]/I");
  tree->Branch("pho_conv_pairinvmass",&pho_conv_pairinvmass,"pho_conv_pairinvmass[pho_n]/F");
  tree->Branch("pho_conv_paircotthetasep",&pho_conv_paircotthetasep,"pho_conv_paircotthetasep[pho_n]/F");
  tree->Branch("pho_conv_eoverp",&pho_conv_eoverp,"pho_conv_eoverp[pho_n]/F");
  tree->Branch("pho_conv_zofprimvtxfromtrks",&pho_conv_zofprimvtxfromtrks,"pho_conv_zofprimvtxfromtrks[pho_n]/F");
  tree->Branch("pho_conv_distofminapproach",&pho_conv_distofminapproach,"pho_conv_distofminapproach[pho_n]/F");
  tree->Branch("pho_conv_dphitrksatvtx",&pho_conv_dphitrksatvtx,"pho_conv_dphitrksatvtx[pho_n]/F");
  tree->Branch("pho_conv_dphitrksatecal",&pho_conv_dphitrksatecal,"pho_conv_dphitrksatecal[pho_n]/F");
  tree->Branch("pho_conv_detatrksatecal",&pho_conv_detatrksatecal,"pho_conv_detatrksatecal[pho_n]/F");
  tree->Branch("pho_conv_tk1_d0",&pho_conv_tk1_d0,"pho_conv_tk1_d0[pho_n]/F");
  tree->Branch("pho_conv_tk1_pout",&pho_conv_tk1_pout,"pho_conv_tk1_pout[pho_n]/F");
  tree->Branch("pho_conv_tk1_pin",&pho_conv_tk1_pin,"pho_conv_tk1_pin[pho_n]/F");
  tree->Branch("pho_conv_tk2_d0",&pho_conv_tk2_d0,"pho_conv_tk2_d0[pho_n]/F");
  tree->Branch("pho_conv_tk2_pout",&pho_conv_tk2_pout,"pho_conv_tk2_pout[pho_n]/F");
  tree->Branch("pho_conv_tk2_pin",&pho_conv_tk2_pin,"pho_conv_tk2_pin[pho_n]/F");
  
  //added by marco
  tree->Branch("pho_conv_tk1_dz",&pho_conv_tk1_dz,"pho_conv_tk1_dz[pho_n]/F");
  tree->Branch("pho_conv_tk2_dz",&pho_conv_tk2_dz,"pho_conv_tk2_dz[pho_n]/F");
  tree->Branch("pho_conv_tk1_dzerr",&pho_conv_tk1_dzerr,"pho_conv_tk1_dzerr[pho_n]/F");
  tree->Branch("pho_conv_tk2_dzerr",&pho_conv_tk2_dzerr,"pho_conv_tk2_dzerr[pho_n]/F");
  tree->Branch("pho_conv_tk1_nh",&pho_conv_tk1_nh,"pho_conv_tk1_nh[pho_n]/I");
  tree->Branch("pho_conv_tk2_nh",&pho_conv_tk2_nh,"pho_conv_tk2_nh[pho_n]/I");
  tree->Branch("pho_conv_chi2",&pho_conv_chi2,"pho_conv_chi2[pho_n]/F");
  tree->Branch("pho_conv_chi2_probability",&pho_conv_chi2_probability,"pho_conv_chi2_probability[pho_n]/F");
  tree->Branch("pho_conv_ch1ch2",&pho_conv_ch1ch2,"pho_conv_ch1ch2[pho_n]/I");
  tree->Branch("pho_conv_validvtx",&pho_conv_validvtx,"pho_conv_validvtx[pho_n]/I");
  tree->Branch("pho_conv_MVALikelihood",&pho_conv_MVALikelihood,"pho_conv_MVALikelihood[pho_n]/I");

  // added by pasquale
  tree->Branch("pho_sipip",&pho_sipip,"pho_sipip[pho_n]/F");
  tree->Branch("pho_sieip",&pho_sieip,"pho_sieip[pho_n]/F");
  tree->Branch("pho_zernike20",&pho_zernike20,"pho_zernike20[pho_n]/F");
  tree->Branch("pho_zernike42",&pho_zernike42,"pho_zernike42[pho_n]/F");
  tree->Branch("pho_e2nd",&pho_e2nd,"pho_e2nd[pho_n]/F");
  tree->Branch("pho_e5x5",&pho_e5x5,"pho_e5x5[pho_n]/F");

  tree->Branch("pho_e2x5right",&pho_e2x5right,"pho_e2x5right[pho_n]/F");
  tree->Branch("pho_e2x5left",&pho_e2x5left,"pho_e2x5left[pho_n]/F");
  tree->Branch("pho_e2x5Top",&pho_e2x5Top,"pho_e2x5Top[pho_n]/F");
  tree->Branch("pho_e2x5bottom",&pho_e2x5bottom,"pho_e2x5bottom[pho_n]/F");
  tree->Branch("pho_eright",&pho_eright,"pho_eright[pho_n]/F");
  tree->Branch("pho_eleft",&pho_eleft,"pho_eleft[pho_n]/F");
  tree->Branch("pho_etop",&pho_etop,"pho_etop[pho_n]/F");
  tree->Branch("pho_ebottom",&pho_ebottom,"pho_ebottom[pho_n]/F");

  //tree->Branch("pho_e2overe9",&pho_e2overe9,"pho_e2overe9[pho_n]/F");
  tree->Branch("pho_seed_severity",&pho_seed_severity,"pho_seed_severity[pho_n]/F");
  tree->Branch("pho_seed_time",&pho_seed_time,"pho_seed_time[pho_n]/F");
  tree->Branch("pho_seed_outoftimechi2",&pho_seed_outoftimechi2,"pho_seed_outoftimechi2[pho_n]/F");
  tree->Branch("pho_seed_chi2",&pho_seed_chi2,"pho_seed_chi2[pho_n]/F");
  tree->Branch("pho_seed_recoflag",&pho_seed_recoflag,"pho_seed_recoflag[pho_n]/F");

  tree->Branch("pho_isconv", &pho_isconv, "pho_isconv[pho_n]/I");
  tree->Branch("pho_residCorrEnergy", &pho_residCorrEnergy, "pho_residCorrEnergy[pho_n]/F");
  tree->Branch("pho_residCorrResn", &pho_residCorrResn, "pho_residCorrResn[pho_n]/F");
  tree->Branch("pho_id", &pho_id, "pho_id[pho_n]/I");
  
  pho_conv_vtx = new TClonesArray("TVector3", MAX_PHOTONS);
  tree->Branch("pho_conv_vtx", "TClonesArray", &pho_conv_vtx, 32000, 0);
  pho_conv_pair_momentum = new TClonesArray("TVector3", MAX_PHOTONS);
  tree->Branch("pho_conv_pair_momentum", "TClonesArray", &pho_conv_pair_momentum, 32000, 0);
  pho_conv_refitted_momentum = new TClonesArray("TVector3", MAX_PHOTONS);
  tree->Branch("pho_conv_refitted_momentum", "TClonesArray", &pho_conv_refitted_momentum, 32000, 0);
  pho_conv_vertexcorrected_p4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
  tree->Branch("pho_conv_vertexcorrected_p4", "TClonesArray", &pho_conv_vertexcorrected_p4, 32000, 0);
}

bool GlobePhotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  if (debug_level > 9) 
    std::cout << "GlobePhotons: Start analyze" << std::endl;

  PhotonFixCMS::initialise(iSetup, "4_2");

  // get collections
  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);

  ///////////////
  // take the pi0 rejection info from RECO
  edm::Handle<reco::PhotonPi0DiscriminatorAssociationMap>  map;
  reco::PhotonPi0DiscriminatorAssociationMap::const_iterator mapIter;
  if (!doAodSim) 
    iEvent.getByLabel("piZeroDiscriminators","PhotonPi0DiscriminatorAssociationMap",  map);

  edm::Handle<reco::PhotonCollection> R_PhotonHandle;
  iEvent.getByLabel(photonCollStd, R_PhotonHandle);
  const reco::PhotonCollection R_photons = *(R_PhotonHandle.product());   

  // Transform Track into TransientTrack (needed by the Vertex fitter)
  edm::ESHandle<TransientTrackBuilder> theTTkBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTkBuilder);

  //// Get the Out In CKF tracks from conversions 
  bool validTrackInputs=true;
  edm::Handle<reco::TrackCollection> outInTrkHandle;
  iEvent.getByLabel("ckfOutInTracksFromConversions",  outInTrkHandle);

  if (!outInTrkHandle.isValid()) {
    std::cout << "Error! Can't get the conversionOITrack " << "\n";
    validTrackInputs=false;
    if (debug_level > 9)
      std::cout  << "ConvertedPhotonProducer  outInTrack collection size " << (*outInTrkHandle).size() << "\n";
  }

  //// Get the association map between CKF Out In tracks and the SC where they originated
  edm::Handle<reco::TrackCaloClusterPtrAssociation> outInTrkSCAssocHandle;
  iEvent.getByLabel( "ckfOutInTracksFromConversions" , "outInTrackSCAssociationCollection", outInTrkSCAssocHandle);
  if (!outInTrkSCAssocHandle.isValid()) {
    //  std::cout << "Error! Can't get the product " <<  outInTrackSCAssociationCollection_.c_str() <<"\n";
    validTrackInputs=false;
  }

  std::vector<reco::TransientTrack> t_outInTrk = ( *theTTkBuilder).build(outInTrkHandle );

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamSpotColl, bsHandle);
  const reco::BeamSpot &thebs = *bsHandle.product();
  
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(convertedPhotonColl, hConversions);
  
  //edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel(electronColl, hElectrons);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoCollection, rhoHandle);
  rho = *(rhoHandle.product());
  
  //edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByLabel(vtxCollection, hVertex);
  iEvent.getByLabel(tkCollection, tkHandle);

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;

  // FOR PF ISOLATION
  edm::Handle<reco::PFCandidateCollection> pfCollection;
  iEvent.getByLabel(pfColl, pfCollection);

  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel("pfSelectedPhotons", pfHandle);
 
  std::vector< edm::Handle< edm::ValueMap<double> > > isolationValues03(inputTagIsoVals03_.size());
  for (size_t j = 0; j<inputTagIsoVals03_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoVals03_[j], isolationValues03[j]);
  }
  
  std::vector< edm::Handle< edm::ValueMap<double> > > isolationValues04(inputTagIsoVals04_.size());
  for (size_t j = 0; j<inputTagIsoVals04_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoVals04_[j], isolationValues04[j]);
  }

  if (debug_level > 9) {
    std::cout << "GlobePhotons: Start analyze" << std::endl;
  }

  // Initialise the Correction Scheme
  fEtaCorr->init(iSetup);
  CrackCorr->init(iSetup);
  LocalCorr->init(iSetup);

  edm::Handle<reco::SuperClusterCollection> superClustersHybridH; 
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 
  edm::Handle<reco::BasicClusterShapeAssociationCollection> hybridClusterShapeBarrelH; 
  edm::Handle<reco::BasicClusterShapeAssociationCollection> basicClusterShapeEndcapH; 

  iEvent.getByLabel(barrelHybridClusterShapeColl, hybridClusterShapeBarrelH);
  iEvent.getByLabel(endcapBasicClusterShapeColl, basicClusterShapeEndcapH);

  iEvent.getByLabel(hybridSuperClusterColl,superClustersHybridH);
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  if (debug_level > 9) {
    std::cout << "GlobePhotons: Photon collection size: "<< phoH->size() << std::endl;

    std::cout << "GlobePhotons: superClustersEndcap collection size: "<< superClustersEndcapH->size() << std::endl;
  }

  // now have collections
  pho_p4->Clear();
  pho_calopos->Clear();
  pho_conv_vtx->Clear();
  pho_conv_pair_momentum->Clear();
  pho_conv_refitted_momentum->Clear();
  pho_conv_vertexcorrected_p4->Clear();
  pho_pfiso_mycharged03->clear();
  pho_pfiso_mycharged04->clear();
  pho_pfiso_mycharged03_noveto->clear();
  pho_pfiso_mycharged04_noveto->clear();

  pho_n = 0;

  if(debug_level>9)std::cout << "GlobePhotons: photons" << std::endl;

  for(unsigned int iPho = 0; iPho < phoH->size(); ++iPho) {
    if (pho_n >= MAX_PHOTONS) {
      std::cout << "GlobePhotons: WARNING TOO MANY PHOTONS: " << phoH->size() << " (allowed " << MAX_PHOTONS << ")" << std::endl;
      break;
    }
    
    reco::PhotonRef localPho(phoH, iPho);

    if(gCUT->cut(*localPho)) 
      continue;

    new ((*pho_p4)[pho_n]) TLorentzVector();
    new ((*pho_calopos)[pho_n]) TVector3();
    new ((*pho_conv_vtx)[pho_n]) TVector3();
    new ((*pho_conv_pair_momentum)[pho_n]) TVector3();
    new ((*pho_conv_refitted_momentum)[pho_n]) TVector3();
    new ((*pho_conv_vertexcorrected_p4)[pho_n]) TLorentzVector();

    if(debug_level>9)
      std::cout << "GlobePhotons: -21 "<< std::endl;

    ((TLorentzVector *)pho_p4->At(pho_n))->SetXYZT(localPho->px(), localPho->py(), localPho->pz(), localPho->energy());
    ((TVector3 *)pho_calopos->At(pho_n))->SetXYZ(localPho->caloPosition().x(), localPho->caloPosition().y(), localPho->caloPosition().z());

    reco::SuperClusterRef theClus=localPho->superCluster();

    pho_h[pho_n] = hoeCalculator(&(*(theClus->seed())), geometry, iEvent, iSetup);
    pho_hoe[pho_n] = localPho->hadronicOverEm();
    pho_scind[pho_n] = -1;

    // PHOTON ID
    pho_id[pho_n] = PhotonID(localPho, 4, reco::VertexRef(hVertex, 0), false);

    // Residual corrections
    PhotonFixCMS ResidCorrector(*localPho);
    pho_residCorrEnergy[pho_n] = ResidCorrector.fixedEnergy();
    pho_residCorrResn[pho_n] = ResidCorrector.sigmaEnergy();

    int index = 0;

    for(int isuperClusterType=0; isuperClusterType<3; ++isuperClusterType) {
      if (isuperClusterType == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j){

          reco::SuperClusterRef sc(superClustersHybridH, j);

          //apply cuts
          if(gCUT->cut(*sc))continue;
          //passed cuts
          if (&(*localPho->superCluster()) == &(*sc)) {

            pho_scind[pho_n] = index;
            break;
          }
          index++;
        }
      }

      if (isuperClusterType == 2) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j){

          reco::SuperClusterRef sc(superClustersEndcapH, j);
          //apply cuts
          if(gCUT->cut(*sc))continue;
          //passed cuts

          if (&(*(localPho->superCluster())) == &(*sc)) {
            pho_scind[pho_n] = index;
            break;
          }
          index++;
        }
      }
    }

    DetId id=localPho->superCluster()->seed()->hitsAndFractions()[0].first;

    bool isBarrel=(id.subdetId() == EcalBarrel);
    pho_barrel[pho_n]=(Int_t)isBarrel;


    // Correction Schemes,
    for (int corr_iter=1;corr_iter<=5;corr_iter++) 
      pho_feta[pho_n][corr_iter] = fEtaCorr->getValue(*(localPho->superCluster()),corr_iter);
     
    pho_crackcorr[pho_n] = CrackCorr->getValue(*(localPho->superCluster()));
    if (isBarrel) pho_localcorr[pho_n] = LocalCorr->getValue(*(localPho->superCluster()),0);
    else pho_localcorr[pho_n] = 1.;

    // Rech-Hits related
    EcalClusterLazyTools lazyTool(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl );   
    edm::Handle<EcalRecHitCollection> prechits;
    iEvent.getByLabel( (localPho->isEB() ? ecalHitEBColl : ecalHitEEColl) ,prechits );
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);

    const reco::CaloClusterPtr  seed_clu = localPho->superCluster()->seed();
    EcalRecHitCollection::const_iterator seedcry_rh = prechits->find( id );
    
    //fiducial flags
    pho_isEB[pho_n] = localPho->isEB();
    pho_isEE[pho_n] = localPho->isEE();
    pho_isEBGap[pho_n] = localPho->isEBGap();
    pho_isEEGap[pho_n] = localPho->isEEGap();
    pho_isEBEEGap[pho_n] = localPho->isEBEEGap();


    //shower shape variables
    pho_see[pho_n] = localPho->sigmaEtaEta();
    pho_sieie[pho_n] = localPho->sigmaIetaIeta();
    pho_e1x5[pho_n] = localPho->e1x5();
    pho_e2x5[pho_n] = localPho->e2x5();
    pho_e3x3[pho_n] = localPho->e3x3();
    pho_e5x5[pho_n] = localPho->e5x5();
    pho_emaxxtal[pho_n] = localPho->maxEnergyXtal();
    pho_hoe[pho_n] = localPho->hadronicOverEm();
    pho_h1oe[pho_n] = localPho->hadronicDepth1OverEm();
    pho_h2oe[pho_n] = localPho->hadronicDepth2OverEm();
    pho_r1x5[pho_n] = localPho->r1x5();
    pho_r2x5[pho_n] = localPho->r2x5();
    pho_r9[pho_n] = localPho->r9();

    // Added by Aris - Begin
    if (!doAodSim) {
      int R_nphot = 0;
      float nn = -1.;
      pho_pi0disc[pho_n] = nn;
      for( reco::PhotonCollection::const_iterator  R_phot_iter = R_photons.begin(); R_phot_iter != R_photons.end(); R_phot_iter++) { 
        mapIter = map->find(edm::Ref<reco::PhotonCollection>(R_PhotonHandle,R_nphot));
        if(mapIter!=map->end()) {
          nn = mapIter->val;
        }
        if(localPho->p4() == R_phot_iter->p4()) pho_pi0disc[pho_n] = nn;
        R_nphot++;              
      }
      
      int iTrk=0;
      bool ConvMatch = false;
      for( std::vector<reco::TransientTrack>::iterator  iTk =  t_outInTrk.begin(); iTk !=  t_outInTrk.end(); iTk++) {
        edm::Ref<reco::TrackCollection> trackRef(outInTrkHandle, iTrk );    
        iTrk++;
	
        const reco::CaloClusterPtr  aClus = (*outInTrkSCAssocHandle)[trackRef];
	
	float conv_SC_et = aClus->energy()/cosh(aClus->eta());
	float conv_SC_eta = aClus->eta(); float conv_SC_phi = aClus->phi(); 
	
	if(localPho->superCluster()->position() == aClus->position()) {
	  ConvMatch = true;	      
	  if (debug_level > 9) {
	    std::cout <<  " ---> ConversionTrackPairFinder track from handle hits " 
		      << trackRef->recHitsSize() << " inner pt  " 
		      << sqrt(iTk->track().innerMomentum().perp2()) << "\n";  
	    std::cout << " ---> ConversionTrackPairFinder  Out In track belonging to SC with (Et,eta,phi) (" 
		      << conv_SC_et << "," << conv_SC_eta << "," <<  conv_SC_phi << ")\n"; 
	  }
	} 
      }
      
      pho_IsConvOutIn[pho_n] = (int)ConvMatch; 
    }
    // Added by Aris - End


    // more cluster shapes from Lazy Tools
    std::vector<float> viCov;
    viCov = lazyTool.localCovariances(*seed_clu);
    pho_sipip[pho_n] = viCov[2];
    pho_sieip[pho_n] = viCov[1];
    pho_zernike20[pho_n] = lazyTool.zernike20(*seed_clu);
    pho_zernike42[pho_n] = lazyTool.zernike42(*seed_clu);
    pho_e2nd[pho_n] = lazyTool.e2nd(*seed_clu);
    pho_e2x5right[pho_n] = lazyTool.e2x5Right(*seed_clu);
    pho_e2x5left[pho_n] = lazyTool.e2x5Left(*seed_clu);
    pho_e2x5Top[pho_n] = lazyTool.e2x5Top(*seed_clu);
    pho_e2x5bottom[pho_n] = lazyTool.e2x5Bottom(*seed_clu);
    pho_eright[pho_n] = lazyTool.eRight(*seed_clu);
    pho_eleft[pho_n] = lazyTool.eLeft(*seed_clu);
    pho_etop[pho_n] = lazyTool.eTop(*seed_clu);
    pho_ebottom[pho_n] = lazyTool.eBottom(*seed_clu);

    //spike-ID
    //pho_e2overe9[pho_n] = EcalSeverityLevelAlgo::E2overE9( id, *prechits, 5.0, 0.0);
    edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
    const EcalSeverityLevelAlgo* sevLevel = sevlv.product();
    pho_seed_severity[pho_n] = sevLevel->severityLevel(id, *prechits);

    pho_seed_time[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->time() : 999.;
    pho_seed_outoftimechi2[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->outOfTimeChi2() : 999.;
    pho_seed_chi2[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->chi2() : 999.;
    pho_seed_recoflag[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->recoFlag() : 999.;

    //isolation variables
    int myIndex = -1;

    float dRmin = 0.2;
    for (unsigned int j=0; j<pfHandle->size(); ++j) {

      reco::PFCandidatePtr temp(pfHandle, j);
      float dR = deltaR(localPho->superCluster()->eta(), localPho->superCluster()->phi(), temp->superClusterRef()->eta(), temp->superClusterRef()->phi());
      if (dR < dRmin) 
	myIndex = j;
    }

    if (myIndex != -1) {
      reco::PFCandidatePtr pfCandidate(pfHandle, myIndex);
	
      pho_pfiso_photon03[pho_n]  = (*isolationValues03[1])[pfCandidate];
      pho_pfiso_neutral03[pho_n] = (*isolationValues03[2])[pfCandidate];
      pho_pfiso_photon03_noveto[pho_n]  = (*isolationValues03[3])[pfCandidate];

      pho_pfiso_photon04[pho_n]  = (*isolationValues04[1])[pfCandidate];
      pho_pfiso_neutral04[pho_n] = (*isolationValues04[2])[pfCandidate];
      pho_pfiso_photon04_noveto[pho_n]  = (*isolationValues04[3])[pfCandidate];
      
      pho_pfiso_mycharged03->push_back(pfTkIsoWithVertex(pfCandidate, pfCollection.product(), 0.3, 0.02)); 
      pho_pfiso_mycharged04->push_back(pfTkIsoWithVertex(pfCandidate, pfCollection.product(), 0.4, 0.02)); 

      pho_pfiso_mycharged03_noveto->push_back(pfTkIsoWithVertex(pfCandidate, pfCollection.product(), 0.3, 0)); 
      pho_pfiso_mycharged04_noveto->push_back(pfTkIsoWithVertex(pfCandidate, pfCollection.product(), 0.4, 0)); 
				       
    } else {
      math::XYZVector vCand(localPho->caloPosition().x(), localPho->caloPosition().y(), localPho->caloPosition().z());
      pho_pfiso_photon03[pho_n]  = pfEcalIso(vCand, pfCollection.product(), 0.3, 0.045, 0.00, 0.0, 0.08, 0.1);
      pho_pfiso_neutral03[pho_n] = pfHcalIso(vCand, pfCollection.product(), 0.3, 0.00);
      pho_pfiso_photon03_noveto[pho_n]  = pfEcalIso(vCand, pfCollection.product(), 0.4, 0.045, 0.0, 0.0, 0.08, 0.1);
      
      pho_pfiso_photon04[pho_n]  = pfEcalIso(vCand, pfCollection.product(), 0.4, 0.045, 0.00, 0.0, 0.08, 0.1);
      pho_pfiso_neutral04[pho_n] = pfHcalIso(vCand, pfCollection.product(), 0.4, 0.00);
      pho_pfiso_photon04_noveto[pho_n] = pfEcalIso(vCand, pfCollection.product(), 0.4, 0.045, 0.0, 0.0, 0.08, 0.1);

      pho_pfiso_mycharged03->push_back(pfTkIsoWithVertex(vCand, pfCollection.product(), 0.3, 0.02)); 
      pho_pfiso_mycharged04->push_back(pfTkIsoWithVertex(vCand, pfCollection.product(), 0.4, 0.02)); 

      pho_pfiso_mycharged03_noveto->push_back(pfTkIsoWithVertex(vCand, pfCollection.product(), 0.3, 0.02)); 
      pho_pfiso_mycharged04_noveto->push_back(pfTkIsoWithVertex(vCand, pfCollection.product(), 0.4, 0.02)); 
    }

    pho_ecalsumetconedr04[pho_n] = localPho->ecalRecHitSumEtConeDR04();
    pho_hcalsumetconedr04[pho_n] = localPho->hcalTowerSumEtConeDR04();
    pho_hcal1sumetconedr04[pho_n] = localPho->hcalDepth1TowerSumEtConeDR04();
    pho_hcal2sumetconedr04[pho_n] = localPho->hcalDepth2TowerSumEtConeDR04();
    pho_trksumptsolidconedr04[pho_n] = localPho->trkSumPtSolidConeDR04();
    pho_trksumpthollowconedr04[pho_n] = localPho->trkSumPtHollowConeDR04();
    pho_ntrksolidconedr04[pho_n] = localPho->nTrkSolidConeDR04();
    pho_ntrkhollowconedr04[pho_n] = localPho->nTrkHollowConeDR04();
    pho_ecalsumetconedr03[pho_n] = localPho->ecalRecHitSumEtConeDR03();
    pho_hcalsumetconedr03[pho_n] = localPho->hcalTowerSumEtConeDR03();
    pho_hcal1sumetconedr03[pho_n] = localPho->hcalDepth1TowerSumEtConeDR03();
    pho_hcal2sumetconedr03[pho_n] = localPho->hcalDepth2TowerSumEtConeDR03();
    pho_trksumptsolidconedr03[pho_n] = localPho->trkSumPtSolidConeDR03();
    pho_trksumpthollowconedr03[pho_n] = localPho->trkSumPtHollowConeDR03();
    pho_ntrksolidconedr03[pho_n] = localPho->nTrkSolidConeDR03();
    pho_ntrkhollowconedr03[pho_n] = localPho->nTrkHollowConeDR03();

    bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(localPho->superCluster(), hElectrons, hConversions, thebs.position());
    pho_isconv[pho_n] = int(passelectronveto);

    //other variables
    pho_haspixseed[pho_n] = localPho->hasPixelSeed();

    pho_hasconvtks[pho_n] = localPho->hasConversionTracks();
    pho_nconv[pho_n] = localPho->conversions().size();

    pho_conv_ntracks[pho_n]=0;
    pho_conv_pairinvmass[pho_n]=-999.;
    pho_conv_paircotthetasep[pho_n]=-999.;
    pho_conv_eoverp[pho_n]=-999.;
    pho_conv_zofprimvtxfromtrks[pho_n]=-999.;
    pho_conv_distofminapproach[pho_n]=-999.;
    pho_conv_dphitrksatvtx[pho_n]=-999.;
    pho_conv_dphitrksatecal[pho_n]=-999.;
    pho_conv_detatrksatecal[pho_n]=-999.;
    pho_conv_tk1_d0[pho_n]=-999.;
    pho_conv_tk1_pout[pho_n]=-999.;
    pho_conv_tk1_pin[pho_n]=-999.;
    pho_conv_validvtx[pho_n]=0;

    ((TVector3 *)pho_conv_vtx->At(pho_n))->SetXYZ(-999, -999, -999);
    ((TVector3 *)pho_conv_pair_momentum->At(pho_n))->SetXYZ(-999, -999, -999);
    ((TVector3 *)pho_conv_refitted_momentum->At(pho_n))->SetXYZ(-999, -999, -999);

    if (debug_level>9) std::cout << "Looking For Valid Conversion" << std::endl;
    if (localPho->hasConversionTracks()) {
      if (debug_level>9) std::cout << "Has Conversion Tracks" << std::endl;

      //reco::ConversionRef conv(localPho->conversions(),0);

      reco::ConversionRefVector conversions = localPho->conversions();

      reco::ConversionRef conv = conversions[0];

      if (debug_level>9) std::cout << "Checking Vertex Validity" << std::endl;
      pho_conv_validvtx[pho_n]=conv->conversionVertex().isValid();

      if (!pho_conv_validvtx[pho_n]) {
        if (debug_level>9) std::cout << "Invalid Conversion" << std::endl;
        ((TLorentzVector *)pho_conv_vertexcorrected_p4->At(pho_n))->SetXYZT(localPho->px(), localPho->py(), localPho->pz(), localPho->energy());
        pho_n++;
        continue;
      }
      
      if (debug_level>9) std::cout << "Starting Conversion Loop" << std::endl;

      for (unsigned int i=0; i<conversions.size(); i++) {
	
        //	reco::ConversionRef 
        conv=conversions[i];

        reco::Vertex vtx=conv->conversionVertex();

        //std::cout<<"Marco Conv vtx: R, x, y, z "<<vtx.position().R()<<" "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<" "<<vtx.chi2()<<std::endl;
        ((TVector3 *)pho_conv_vtx->At(pho_n))->SetXYZ(vtx.x(), vtx.y(), vtx.z());
        ((TVector3 *)pho_conv_pair_momentum->At(pho_n))->SetXYZ(conv->pairMomentum().x(), conv->pairMomentum().y(), conv->pairMomentum().z());
        ((TVector3 *)pho_conv_refitted_momentum->At(pho_n))->SetXYZ(conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z());
        
        pho_conv_chi2[pho_n]=vtx.chi2();
        pho_conv_chi2_probability[pho_n]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
        pho_conv_ntracks[pho_n]=conv->nTracks();
        pho_conv_MVALikelihood[pho_n]=conv->MVAout();

        if(pho_conv_ntracks[pho_n]) {
	  const std::vector<edm::RefToBase<reco::Track> > tracks = conv->tracks();
	  for (unsigned int i=0; i<tracks.size(); i++) {
            if(i==0) {
              pho_conv_tk1_dz[pho_n]=tracks[i]->dz();
              pho_conv_tk1_dzerr[pho_n]=tracks[i]->dzError();
              pho_conv_ch1ch2[pho_n]=tracks[i]->charge();
            } else if(i==1) {
              pho_conv_tk2_dz[pho_n]=tracks[i]->dz();
              pho_conv_tk2_dzerr[pho_n]=tracks[i]->dzError();
              pho_conv_ch1ch2[pho_n]*=tracks[i]->charge();
            }
          }
        }
      }

      pho_conv_pairinvmass[pho_n]=conv->pairInvariantMass();
      pho_conv_paircotthetasep[pho_n]=conv->pairCotThetaSeparation();
      pho_conv_eoverp[pho_n]=conv->EoverPrefittedTracks();
      pho_conv_zofprimvtxfromtrks[pho_n]=conv->zOfPrimaryVertexFromTracks();
      pho_conv_distofminapproach[pho_n]=conv->distOfMinimumApproach();
      pho_conv_dphitrksatvtx[pho_n]=conv->dPhiTracksAtVtx();
      pho_conv_dphitrksatecal[pho_n]=conv->dPhiTracksAtEcal();
      pho_conv_detatrksatecal[pho_n]=conv->dEtaTracksAtEcal();

      if(conv->tracks().size() > 0) {
        pho_conv_tk1_d0[pho_n]=conv->tracksSigned_d0()[0];
        pho_conv_tk1_pout[pho_n]=sqrt(conv->tracksPout()[0].Mag2());
        pho_conv_tk1_pin[pho_n]=sqrt(conv->tracksPin()[0].Mag2());
      }

      if(conv->tracks().size() > 1) {
        pho_conv_tk2_d0[pho_n]=conv->tracksSigned_d0()[1];
        pho_conv_tk2_pout[pho_n]=sqrt(conv->tracksPout()[1].Mag2());
        pho_conv_tk2_pin[pho_n]=sqrt(conv->tracksPin()[1].Mag2());
      }

      reco::Photon temp = *localPho;
      temp.setVertex(math::XYZPoint(0.0, 0.0, pho_conv_zofprimvtxfromtrks[pho_n]));
      ((TLorentzVector *)pho_conv_vertexcorrected_p4->At(pho_n))->SetXYZT(temp.px(), temp.py(), temp.pz(), temp.energy());
    }

    pho_n++;
  }

  if(debug_level>9)
  std::cout << "End Photon" << std::endl;

  return true;
}

float GlobePhotons::hoeCalculator(const reco::BasicCluster* clus, const CaloGeometry& geometry,
				  const edm::Event& e , const edm::EventSetup& c) {
  
  float h = 0.;

  GlobalPoint pclu(clus->x(),clus->y(),clus->z());

  edm::Handle< HBHERecHitCollection > hbhe ;
  e.getByLabel(hcalHitColl, hbhe);
  const HBHERecHitCollection* hithbhe_ = hbhe.product();
 
  const CaloSubdetectorGeometry *geometry_p ; 
  geometry_p = geometry.getSubdetectorGeometry (DetId::Hcal,4) ;

  DetId hcalDetId ;
  hcalDetId = geometry_p->getClosestCell(pclu) ;

  CaloRecHitMetaCollection f;
  f.add(hithbhe_);
  CaloRecHitMetaCollection::const_iterator iterRecHit; 
  iterRecHit = f.find(hcalDetId) ;
  if (iterRecHit!=f.end()) {
    h = iterRecHit->energy() ;
  }
 
  return h;
}

LorentzVector GlobePhotons::VertexCorrectedP4Hgg(reco::PhotonRef photon, reco::VertexRef vtx) {
  
  if (vtx.isNull())
    return photon->p4();
  
  math::XYZPointF ecalPosition = photon->caloPosition();
  math::XYZPoint vertexPosition = vtx->position();
  math::XYZVector u ((ecalPosition - vertexPosition).x(), (ecalPosition - vertexPosition).y(), (ecalPosition - vertexPosition).z());

  double energy = photon->energy();  
  u = u.Unit();
  LorentzVector vertexCorrectedP4(energy*u.x(), energy*u.y(), energy*u.z(), energy);

  return vertexCorrectedP4;
}

int GlobePhotons::PhotonID(reco::PhotonRef photon, int ncat, reco::VertexRef vtx, bool doSublead, int diphoind)  {

  int cutlevelpassed = -1;

  int n_r9_categories = -1;
  int n_eta_categories = -1;
  if(ncat == 6) {
    n_r9_categories = 3;
    n_eta_categories = 2;
  } else if(ncat == 4) {
    n_r9_categories = 2;
    n_eta_categories = 2;
  }

  LorentzVector pholv = VertexCorrectedP4Hgg(photon, vtx);

  float val_sieie = photon->sigmaIetaIeta();
  float val_hoe = photon->hadronicOverEm();
  float val_r9 = photon->r9();
  float val_pixel = (float)photon->hasPixelSeed();
  float	val_drtotk = DeltaRToTrackHgg(photon, hElectrons, 0);
  float	val_tkiso =  SumTrackPtInConeHgg(photon, vtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
  float val_tkisobad = WorstSumTrackPtInConeHgg(photon, 0., 0.40, 0.02, 0.0, 1.0, 0.1);

  if(!vtx.isNull()) {
    if(diphoind != -1) {
      if(!doSublead) {
	val_drtotk = DeltaRToTrackHgg(photon, hElectrons, 0);
	val_tkiso =  SumTrackPtInConeHgg(photon, vtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
	val_tkisobad = WorstSumTrackPtInConeHgg(photon, 0., 0.40, 0.02, 0.0, 1.0, 0.1);
      }
      else {
	val_drtotk = DeltaRToTrackHgg(photon, hElectrons, 0);
	val_tkiso = SumTrackPtInConeHgg(photon, vtx, 0, 0.30, 0.02, 0.0, 1.0, 0.1);
	val_tkisobad = WorstSumTrackPtInConeHgg(photon, 0., 0.40, 0.02, 0.0, 1.0, 0.1);
      }
    }
  }

  float rhofacbad = 0.52, rhofac = 0.17;

  float val_isosumoet = (val_tkiso + photon->ecalRecHitSumEtConeDR03() +
			 photon->hcalTowerSumEtConeDR04() - rho*rhofac)*50./pholv.Et();
  float val_isosumoetbad = (val_tkisobad + photon->ecalRecHitSumEtConeDR04() + 
			    photon->hcalTowerSumEtConeDR04() - rho*rhofacbad)*50./pholv.Et();
  float val_trkisooet = (val_tkiso)*50./pholv.Et();
  
  if (debug_level > 9) {
    std::cout << "val_isosumoet "    << val_isosumoet << std::endl;
    std::cout << "val_isosumoetbad " << val_isosumoetbad << std::endl;
    std::cout << "val_trkisooet "    << val_trkisooet << std::endl;
    std::cout << "val_sieie "        << val_sieie << std::endl;
    std::cout << "val_hoe "          << val_hoe << std::endl;
    std::cout << "val_r9 "           << val_r9 << std::endl;
    std::cout << "val_drtotk "       << val_drtotk << std::endl;
    std::cout << "val_pixel "        << val_pixel << std::endl;    
    std::cout << val_tkiso << " " << val_tkisobad << " " << pholv.Et() << std::endl;
  }
  
  float val_eta = fabs(photon->caloPosition().eta());

  if (ncat == 4) { 
    cutlevelpassed = photonCutLevel4(val_r9, val_eta, val_isosumoet, val_isosumoetbad, val_trkisooet, val_sieie, val_hoe, val_drtotk); 
  } else if(ncat == 6) { 
    cutlevelpassed = photonCutLevel6(val_r9, val_eta, val_isosumoet, val_isosumoetbad, val_trkisooet, val_sieie, val_hoe, val_drtotk);
  } else { 
    std::cerr << "Photon selection for " << ncat << " categories does not exist" << std::endl; 
  }
  
  return cutlevelpassed;
}

int GlobePhotons::photonCutLevel4(float r9, float eta, float isosumoet, float isosumoetbad, float tkisosumoet, float sieie, float hoe, float drtotk) {

  int ccat=0;
  if (fabs(eta) > 1.479) 
    ccat = 2;                       //   4 cats
  if (r9 < 0.94) 
    ccat++;

  int lev=0;
  
  while(lev < 11 &&
	isosumoet    <= cutsubleadisosumoet[lev][ccat] &&
	isosumoetbad <= cutsubleadisosumoetbad[lev][ccat] &&
	tkisosumoet  <= cutsubleadtrkisooetom[lev][ccat] &&
	sieie        <= cutsubleadsieie[lev][ccat] &&
	hoe          <= cutsubleadhovere[lev][ccat] &&
	r9           >= cutsubleadr9[lev][ccat] &&
	drtotk       >= cutsublead_drtotk[lev][ccat])
    
    lev++;

  return lev;
}

int GlobePhotons::photonCutLevel6(float r9, float eta, float isosumoet, float isosumoetbad, float tkisosumoet, float sieie, float hoe, float drtotk) {
  
  int ccat = 0;
  if (fabs(eta) > 1.479) 
    ccat = 3;
  if (r9 < 0.94) 
    ccat++;
  if (r9 < 0.90)
    ccat++;

  int lev=0;
  while(lev < 11 && 
	isosumoet    <= cutsubleadisosumoet6c[lev][ccat] &&
    	isosumoetbad <= cutsubleadisosumoetbad6c[lev][ccat] &&
	tkisosumoet  <= cutsubleadtrkisooetom6c[lev][ccat] &&
	sieie        <= cutsubleadsieie6c[lev][ccat] &&
	hoe          <= cutsubleadhovere6c[lev][ccat] &&
	r9           >= cutsubleadr96c[lev][ccat] &&
	drtotk       >= cutsublead_drtotk6c[lev][ccat]) 
    lev++;

  return lev;
}

Float_t GlobePhotons::SumTrackPtInConeHgg(reco::PhotonRef photon, reco::VertexRef vtx, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  if (vtx.isNull())
    return -99;
  
  float SumTrackPt = 0;

  for(unsigned int itk = 0; itk < tkHandle->size(); ++itk) {
    
    reco::TrackRef tk(tkHandle, itk);

    if(tk->pt() < PtMin)
      continue;

    double deltaz = fabs(vtx->z() - tk->vz());

    if(deltaz > dzmax)
      continue;
    
    double dxy = tk->dxy(vtx->position());
    if(fabs(dxy) > dxymax)
      continue;

    double deta = fabs(photon->eta() - tk->eta());
    double dphi = fabs(photon->phi() - tk->phi());
    if(dphi > TMath::Pi())
      dphi = TMath::TwoPi() - dphi;

    double deltaR = sqrt(deta*deta + dphi*dphi);
    
    if(deltaR < OuterConeRadius && deltaR >= InnerConeRadius && deta >= EtaStripHalfWidth) {
      SumTrackPt += tk->pt();
    }
  }

  return SumTrackPt;
}

Float_t GlobePhotons::WorstSumTrackPtInConeHgg(reco::PhotonRef photon, Float_t PtMin, Float_t OuterConeRadius, Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) {

  Int_t worstvtxind = -1;
  Float_t maxisosum = -100;
  
  for(unsigned int ivtx = 0; ivtx < hVertex->size(); ++ivtx) {

    reco::VertexRef vtx(hVertex, ivtx);

    Float_t thisvtxisosum = SumTrackPtInConeHgg(photon, vtx, PtMin, OuterConeRadius, InnerConeRadius, EtaStripHalfWidth, dzmax, dxymax);
    if(thisvtxisosum > maxisosum) {
      maxisosum = thisvtxisosum;
      worstvtxind = ivtx;
    }
  }

  return maxisosum;
}

Float_t GlobePhotons::DeltaRToTrackHgg(reco::PhotonRef photon, edm::Handle<reco::GsfElectronCollection> elHandle, int maxlosthits) {

  int elind = -1;
  float eldr = 99.;

  for(unsigned int iel = 0; iel < elHandle->size(); ++iel) {

    reco::GsfElectronRef el(elHandle, iel);

    if(el->gsfTrack()->trackerExpectedHitsInner().numberOfHits() > maxlosthits)
      continue;

    if(el->superCluster() == photon->superCluster()) {
      elind = iel;
    }
  }

  if(elind >= 0) {
    reco::GsfElectronRef el(elHandle, elind);
    eldr = sqrt(pow(el->deltaEtaSuperClusterTrackAtVtx(), 2) + pow(el->deltaPhiSuperClusterTrackAtVtx(), 2));
  }

  return eldr;
}


std::vector<float> GlobePhotons::pfTkIsoWithVertex(const reco::PFCandidatePtr cand, const reco::PFCandidateCollection* forIsolation, float dRmax, float dRveto) {
  
  std::vector<float> result;
  
  for(unsigned int ivtx=0; ivtx<hVertex->size(); ++ivtx) {
    
    reco::VertexRef vtx(hVertex, ivtx);
    math::XYZVector vCand(cand->superClusterRef()->position().x()
      - vtx->x(), 
      cand->superClusterRef()->position().y()
      - vtx->y(), 
      cand->superClusterRef()->position().z()
      - vtx->z());
    
    float sum = 0;
    for(unsigned i=0; i<forIsolation->size(); i++) {
    
      const reco::PFCandidate& pfc = (*forIsolation)[i];

      if(sameParticle(cand, pfc)) 
      	continue; 
      
      if (pfc.particleId() != reco::PFCandidate::h)
	continue;
      
      if (pfc.pt() < 1.)
	continue;
      
      float dz = fabs(pfc.vz() - vtx->z());
      if (dz > 10)
	continue;
      
      math::XYZVector pvi(pfc.momentum());
      float dR = deltaR(vCand.Eta(), vCand.Phi(), pvi.Eta(), pvi.Phi());
      
      if(dR > dRmax || dR < dRveto)
	continue;
      
      sum += pfc.pt();
    }
    
    result.push_back(sum);
  }
  
  return result;
}

bool GlobePhotons::sameParticle(const reco::PFCandidate& particle1, const reco::PFCandidate& particle2) const {
  
  double smallNumber = 0.02;
  
  if(particle1.particleId() != particle2.particleId()) return false;
  else if( fabs( particle1.energy() - particle2.energy() ) > smallNumber ) return false;
  else if( fabs( particle1.eta() - particle2.eta() ) > smallNumber ) return false;
  else if( fabs( particle1.eta() - particle2.eta() ) > smallNumber ) return false;
  else return true; 
}

std::vector<float> GlobePhotons::pfTkIsoWithVertex(math::XYZVector cand, const reco::PFCandidateCollection* forIsolation, float dRmax, float dRveto) {
  
  std::vector<float> result;
  
  for(unsigned int ivtx=0; ivtx<hVertex->size(); ++ivtx) {
    
    reco::VertexRef vtx(hVertex, ivtx);
    math::XYZVector vCand(cand.x() - vtx->x(), 
			  cand.y() - vtx->y(), 
			  cand.z() - vtx->z());
    
    float sum = 0;
    for(unsigned i=0; i<forIsolation->size(); i++) {
    
      const reco::PFCandidate& pfc = (*forIsolation)[i];
      
      if (pfc.particleId() != reco::PFCandidate::h)
	continue;
      
      if (pfc.pt() < 1.)
	continue;
      
      float dz = fabs(pfc.vz() - vtx->z());
      if (dz > 10)
	continue;
      
      math::XYZVector pvi(pfc.momentum());
      float dR = deltaR(vCand.Eta(), vCand.Phi(), pvi.Eta(), pvi.Phi());
      
      if(dR > dRmax || dR < dRveto)
	continue;
      
      sum += pfc.pt();
    }
    
    result.push_back(sum);
  }
  
  return result;
}


float GlobePhotons::pfEcalIso(math::XYZVector vCand, const reco::PFCandidateCollection* forIsolation, float dRmax, float dRveto, float etaStrip, float phiStrip, float energyBarrel, float energyEndcap) {
  
  float sum = 0;
  for(unsigned i=0; i<forIsolation->size(); i++) {
    
    const reco::PFCandidate& pfc = (*forIsolation)[i];

    if (pfc.particleId() != reco::PFCandidate::gamma)
      continue;
    
    if (fabs(pfc.positionAtECALEntrance().Eta()) < 1.479) {
      dRveto = 0.045;
      if (fabs(pfc.pt()) < energyBarrel)
	continue;
    } else {
      dRveto = 0.070;
      if (fabs(pfc.energy()) < energyEndcap)
	continue;
    }
            
    math::XYZVector pvi(pfc.momentum());
    float dR = deltaR(vCand.Eta(), vCand.Phi(), pvi.Eta(), pvi.Phi());
    float dEta = fabs(vCand.Eta() - pvi.Eta());
    double dPhi = fabs(vCand.Phi() - pvi.Phi());
    if(dPhi > TMath::Pi())
      dPhi = TMath::TwoPi() - dPhi;
		      
    if (dEta < etaStrip)
	continue;

    if (dPhi < phiStrip)
	continue;

    if(dR > dRmax || dR < dRveto)
      continue;
    
    sum += pfc.pt();
  }
  
  return sum;
}

float GlobePhotons::pfHcalIso(math::XYZVector vCand, const reco::PFCandidateCollection* forIsolation, float dRmax, float dRveto) {
  
  float sum = 0;
  for(unsigned i=0; i<forIsolation->size(); i++) {
    
    const reco::PFCandidate& pfc = (*forIsolation)[i];
    
    //if(sameParticle(cand, pfc)) 
    //  continue; 
    
    if (pfc.particleId() != reco::PFCandidate::h0)
      continue;
    
    if (pfc.pt() < 0.5)
      continue;
    
    math::XYZVector pvi(pfc.momentum());
    float dR = deltaR(vCand.Eta(), vCand.Phi(), pvi.Eta(), pvi.Phi());
    
    if(dR > dRmax || dR < dRveto)
      continue;
    
    sum += pfc.pt();
  }
  
  return sum;
}

