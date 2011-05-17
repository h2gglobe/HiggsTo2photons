//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
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

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/EgammaTrackReco/interface/TrackCaloClusterAssociation.h"
/////////////////

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

  // get the Correction Functions
  fEtaCorr  = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrection",iConfig);
  CrackCorr = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection",iConfig);
  LocalCorr = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
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
  
  //tree->Branch("pho_isconv", &pho_isconv, "pho_isconv[pho_n]/I");
  
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
    {
    std::cout << "GlobePhotons: Start analyze" << std::endl;
  }
  // get collections
  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);

  ///////////////
  // take the pi0 rejection info from RECO
  edm::Handle<reco::PhotonPi0DiscriminatorAssociationMap>  map;
  reco::PhotonPi0DiscriminatorAssociationMap::const_iterator mapIter;
  if (!doAodSim) iEvent.getByLabel("piZeroDiscriminators","PhotonPi0DiscriminatorAssociationMap",  map);

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
  }
  if (debug_level > 9)
    std::cout  << "ConvertedPhotonProducer  outInTrack collection size " << (*outInTrkHandle).size() << "\n";

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
  
  edm::Handle<reco::GsfElectronCollection> hElectrons;
  iEvent.getByLabel(electronColl, hElectrons);

  if (debug_level > 9) {
    std::cout  << "ConvertedPhotonProducer  outInTrack association map with SC collection size " << (*outInTrkSCAssocHandle).size() << "\n";
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
  
  pho_n = 0;

  if(debug_level>9)std::cout << "GlobePhotons: photons" << std::endl;

  for( reco::PhotonCollection::const_iterator  iPho = phoH->begin(); iPho != phoH->end(); iPho++) {
    if (pho_n >= MAX_PHOTONS) {
      std::cout << "GlobePhotons: WARNING TOO MANY PHOTONS: " << phoH->size() << " (allowed " << MAX_PHOTONS << ")" << std::endl;
      break;
    }
    reco::Photon localPho = reco::Photon(*iPho);

    if(gCUT->cut(localPho)) continue;

    new ((*pho_p4)[pho_n]) TLorentzVector();
    new ((*pho_calopos)[pho_n]) TVector3();
    new ((*pho_conv_vtx)[pho_n]) TVector3();
    new ((*pho_conv_pair_momentum)[pho_n]) TVector3();
    new ((*pho_conv_refitted_momentum)[pho_n]) TVector3();
    new ((*pho_conv_vertexcorrected_p4)[pho_n]) TLorentzVector();

    if(debug_level>9)std::cout << "GlobePhotons: -21 "<< std::endl;

    ((TLorentzVector *)pho_p4->At(pho_n))->SetXYZT(localPho.px(), localPho.py(), localPho.pz(), localPho.energy());
    ((TVector3 *)pho_calopos->At(pho_n))->SetXYZ(localPho.caloPosition().x(), localPho.caloPosition().y(), localPho.caloPosition().z());

    reco::SuperClusterRef theClus=localPho.superCluster();

    pho_hoe[pho_n]=-1;

    pho_hoe[pho_n] = localPho.hadronicOverEm();
    pho_scind[pho_n] = -1;

    int index = 0;

    for(int isuperClusterType=0; isuperClusterType<3; ++isuperClusterType) {
      if (isuperClusterType == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j){

          reco::SuperClusterRef sc(superClustersHybridH, j);

          //apply cuts
          if(gCUT->cut(*sc))continue;
          //passed cuts
          if (&(*localPho.superCluster()) == &(*sc)) {

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

          if (&(*(localPho.superCluster())) == &(*sc)) {
            pho_scind[pho_n] = index;
            break;
          }
          index++;
        }
      }
    }

    DetId id=localPho.superCluster()->seed()->hitsAndFractions()[0].first;

    bool isBarrel=(id.subdetId() == EcalBarrel);
    pho_barrel[pho_n]=(Int_t)isBarrel;


    // Correction Schemes,
    for (int corr_iter=1;corr_iter<=5;corr_iter++) 
      pho_feta[pho_n][corr_iter] = fEtaCorr->getValue(*(localPho.superCluster()),corr_iter);
     
    pho_crackcorr[pho_n] = CrackCorr->getValue(*(localPho.superCluster()));
    if (isBarrel) pho_localcorr[pho_n] = LocalCorr->getValue(*(localPho.superCluster()),0);
    else pho_localcorr[pho_n] = 1.;

    // Rech-Hits related
    EcalClusterLazyTools lazyTool(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl );   
    edm::Handle<EcalRecHitCollection> prechits;
    iEvent.getByLabel( (localPho.isEB() ? ecalHitEBColl : ecalHitEEColl) ,prechits );
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);

    const reco::CaloClusterPtr  seed_clu = localPho.superCluster()->seed();
    EcalRecHitCollection::const_iterator seedcry_rh = prechits->find( id );
    
    //fiducial flags
    pho_isEB[pho_n] = localPho.isEB();
    pho_isEE[pho_n] = localPho.isEE();
    pho_isEBGap[pho_n] = localPho.isEBGap();
    pho_isEEGap[pho_n] = localPho.isEEGap();
    pho_isEBEEGap[pho_n] = localPho.isEBEEGap();


    //shower shape variables
    pho_see[pho_n] = localPho.sigmaEtaEta();
    pho_sieie[pho_n] = localPho.sigmaIetaIeta();
    pho_e1x5[pho_n] = localPho.e1x5();
    pho_e2x5[pho_n] = localPho.e2x5();
    pho_e3x3[pho_n] = localPho.e3x3();
    pho_e5x5[pho_n] = localPho.e5x5();
    pho_emaxxtal[pho_n] = localPho.maxEnergyXtal();
    pho_hoe[pho_n] = localPho.hadronicOverEm();
    pho_h1oe[pho_n] = localPho.hadronicDepth1OverEm();
    pho_h2oe[pho_n] = localPho.hadronicDepth2OverEm();
    pho_r1x5[pho_n] = localPho.r1x5();
    pho_r2x5[pho_n] = localPho.r2x5();
    pho_r9[pho_n] = localPho.r9();

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
        if(iPho->p4() == R_phot_iter->p4()) pho_pi0disc[pho_n] = nn;
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
	
	if((*iPho).superCluster()->position() == aClus->position()) {
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
    pho_ecalsumetconedr04[pho_n] = localPho.ecalRecHitSumEtConeDR04();
    pho_hcalsumetconedr04[pho_n] = localPho.hcalTowerSumEtConeDR04();
    pho_hcal1sumetconedr04[pho_n] = localPho.hcalDepth1TowerSumEtConeDR04();
    pho_hcal2sumetconedr04[pho_n] = localPho.hcalDepth2TowerSumEtConeDR04();
    pho_trksumptsolidconedr04[pho_n] = localPho.trkSumPtSolidConeDR04();
    pho_trksumpthollowconedr04[pho_n] = localPho.trkSumPtHollowConeDR04();
    pho_ntrksolidconedr04[pho_n] = localPho.nTrkSolidConeDR04();
    pho_ntrkhollowconedr04[pho_n] = localPho.nTrkHollowConeDR04();
    pho_ecalsumetconedr03[pho_n] = localPho.ecalRecHitSumEtConeDR03();
    pho_hcalsumetconedr03[pho_n] = localPho.hcalTowerSumEtConeDR03();
    pho_hcal1sumetconedr03[pho_n] = localPho.hcalDepth1TowerSumEtConeDR03();
    pho_hcal2sumetconedr03[pho_n] = localPho.hcalDepth2TowerSumEtConeDR03();
    pho_trksumptsolidconedr03[pho_n] = localPho.trkSumPtSolidConeDR03();
    pho_trksumpthollowconedr03[pho_n] = localPho.trkSumPtHollowConeDR03();
    pho_ntrksolidconedr03[pho_n] = localPho.nTrkSolidConeDR03();
    pho_ntrkhollowconedr03[pho_n] = localPho.nTrkHollowConeDR03();

    bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(localPho.superCluster(), hElectrons, hConversions, thebs.position());
    pho_isconv[pho_n] = int(passelectronveto);

    //other variables
    pho_haspixseed[pho_n] = localPho.hasPixelSeed();

    pho_hasconvtks[pho_n] = localPho.hasConversionTracks();
    pho_nconv[pho_n] = localPho.conversions().size();

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
    if (localPho.hasConversionTracks()) {
      if (debug_level>9) std::cout << "Has Conversion Tracks" << std::endl;

      //reco::ConversionRef conv(localPho.conversions(),0);

      reco::ConversionRefVector conversions = localPho.conversions();

      reco::ConversionRef conv = conversions[0];

      if (debug_level>9) std::cout << "Checking Vertex Validity" << std::endl;
      pho_conv_validvtx[pho_n]=conv->conversionVertex().isValid();

      if (!pho_conv_validvtx[pho_n]) {
        if (debug_level>9) std::cout << "Invalid Conversion" << std::endl;
        ((TLorentzVector *)pho_conv_vertexcorrected_p4->At(pho_n))->SetXYZT(localPho.px(), localPho.py(), localPho.pz(), localPho.energy());
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
            }
            else if(i==1) {
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

      localPho.setVertex(math::XYZPoint(0.0,0.0,pho_conv_zofprimvtxfromtrks[pho_n]));
      ((TLorentzVector *)pho_conv_vertexcorrected_p4->At(pho_n))->SetXYZT(localPho.px(), localPho.py(), localPho.pz(), localPho.energy());
    }

    pho_n++;
  }


  if(debug_level>9)
  std::cout << "End Photon" << std::endl;

  return true;
}

double GlobePhotons::getHoE(GlobalPoint pclu, double ecalEnergy, const edm::Event& e , const edm::EventSetup& c ) {

  if ( !theCaloGeom_.isValid() )
    c.get<IdealGeometryRecord>().get(theCaloGeom_) ;
  
  //product the geometry
  theCaloGeom_.product() ;

  //Create a CaloRecHitMetaCollection
  edm::Handle< HBHERecHitCollection > hbhe ;
  e.getByLabel(hcalBEColl,hbhe);
  const HBHERecHitCollection* hithbhe_ = hbhe.product();

  double HoE;
  const CaloGeometry& geometry = *theCaloGeom_ ;
  const CaloSubdetectorGeometry *geometry_p ; 

  geometry_p = geometry.getSubdetectorGeometry (DetId::Hcal,4) ;
  DetId hcalDetId ;
  hcalDetId = geometry_p->getClosestCell(pclu) ;
  double hcalEnergy = 0 ;

  CaloRecHitMetaCollection f;
  f.add(hithbhe_);

  CaloRecHitMetaCollection::const_iterator iterRecHit;
  iterRecHit = f.find(hcalDetId) ;

  if(iterRecHit == f.end()) {
    return 0.;
  }

  hcalEnergy = iterRecHit->energy() ;
  HoE = hcalEnergy/ecalEnergy ;

  return HoE ;
}
