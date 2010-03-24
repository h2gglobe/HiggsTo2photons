#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

GlobePhotons::GlobePhotons(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  debug_level = iConfig.getParameter<int>("Debug_Level");
  doFastSim = iConfig.getParameter<bool>("doFastSim");
  doEgammaSummer09Skim = iConfig.getParameter<bool>("doEgammaSummer09Skim");

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

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobePhotons::defineBranch(TTree* tree) {

  pho_p4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
  pho_calopos = new TClonesArray("TVector3", MAX_PHOTONS);

  tree->Branch("pho_n", &pho_n, "pho_n/I");
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
  tree->Branch("pho_conv_ch1ch2",&pho_conv_ch1ch2,"pho_conv_ch1ch2[pho_n]/I");
  tree->Branch("pho_conv_validvtx",&pho_conv_validvtx,"pho_conv_validvtx[pho_n]/I");
  
  pho_conv_vtx = new TClonesArray("TVector3", MAX_PHOTONS);
  tree->Branch("pho_conv_vtx", "TClonesArray", &pho_conv_vtx, 32000, 0);
}

bool GlobePhotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

#ifdef PIZERODISCRIMINATOR
  edm::Handle<reco::PhotonPi0DiscriminatorAssociationMap>  pi0map;
  reco::PhotonPi0DiscriminatorAssociationMap::const_iterator pi0mapIter;
  if(doConvertedPhoton) //CHECK, should be pizero	
    iEvent.getByLabel("piZeroDiscriminators","PhotonPi0DiscriminatorAssociationMap",  pi0map);
#endif
  if (debug_level > 9) 
    {
    std::cout << "GlobePhotons: Start analyze" << std::endl;
  }
  // get collections
  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);

  if (debug_level > 9) {
    std::cout << "GlobePhotons: Start analyze" << std::endl;
  }

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

    if(debug_level>9)std::cout << "GlobePhotons: -21 "<< std::endl;

    ((TLorentzVector *)pho_p4->At(pho_n))->SetXYZT(localPho.px(), localPho.py(), localPho.pz(), localPho.energy());
    ((TVector3 *)pho_calopos->At(pho_n))->SetXYZ(localPho.caloPosition().x(), localPho.caloPosition().y(), localPho.caloPosition().z());

    //std::cout << "Marco GlobePhotons: et, eta, phi "<< localPho.pt()<<" "<<localPho.eta()<<" "<<localPho.phi()<<" "<<std::endl;

    if(debug_level>9)
      std::cout << "Marco GlobePhotons: -22 "<< std::endl;

    reco::SuperClusterRef theClus=localPho.superCluster();

    //GlobalPoint pSc(localPho.superCluster()->position().x(),  localPho.superCluster()->position().y(), localPho.superCluster()->position().z());
    if(debug_level>9)std::cout << "GlobePhotons: -23 "<< std::endl;

    pho_hoe[pho_n]=-1;
    if(debug_level>9)std::cout << "GlobePhotons: -24 "<< std::endl;

    if(debug_level>9)std::cout << "GlobePhotons: -1 "<< std::endl;

    pho_hoe[pho_n] = localPho.hadronicOverEm();
    if(debug_level>9)std::cout << "GlobePhotons: -2 "<< std::endl;
    pho_scind[pho_n] = -1;
    if(debug_level>9)std::cout << "GlobePhotons: -3 "<< std::endl;

    int index = 0;

    if(debug_level>9)std::cout << "GlobePhotons: 0 "<< std::endl;
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

    if(localPho.hasConversionTracks()) {

      //reco::ConversionRef conv(localPho.conversions(),0);

      reco::ConversionRefVector conversions = localPho.conversions();

      reco::ConversionRef conv;

      for (unsigned int i=0; i<conversions.size(); i++) {
	
	//	reco::ConversionRef 
	conv=conversions[i];

	reco::Vertex vtx=conv->conversionVertex();
	//std::cout<<"Marco Conv vtx: R, x, y, z "<<vtx.position().R()<<" "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<" "<<vtx.chi2()<<std::endl;
	((TVector3 *)pho_conv_vtx->At(pho_n))->SetXYZ(vtx.x(), vtx.y(), vtx.z());
	pho_conv_chi2[pho_n]=vtx.chi2();

	pho_conv_validvtx[pho_n]=0;
	if(vtx.isValid()) {
	  pho_conv_validvtx[pho_n]=1;
	}

	pho_conv_ntracks[pho_n]=conv->nTracks();
      
	if(pho_conv_ntracks[pho_n]) {
	  std::vector<reco::TrackRef> tracks = conv->tracks();
	  for (unsigned int i=0; i<tracks.size(); i++) {
	    if(i==0) {
	      pho_conv_tk1_dz[pho_n]=tracks[i]->dz();
	      pho_conv_tk1_dzerr[pho_n]=tracks[i]->dzError();
	      if(!doEgammaSummer09Skim) pho_conv_tk1_nh[pho_n]=tracks[i]->recHitsSize();
	      else pho_conv_tk1_nh[pho_n]= -1;
	      pho_conv_ch1ch2[pho_n]=tracks[i]->charge();
	    }
	    else if(i==1) {
	      pho_conv_tk2_dz[pho_n]=tracks[i]->dz();
	      pho_conv_tk2_dzerr[pho_n]=tracks[i]->dzError();
	      if(!doEgammaSummer09Skim) pho_conv_tk2_nh[pho_n]=tracks[i]->recHitsSize();
	      else pho_conv_tk2_nh[pho_n]= -1;
	      pho_conv_ch1ch2[pho_n]*=tracks[i]->charge();
	    }
	  }
	}
      }
      
      pho_conv_pairinvmass[pho_n]=conv->pairInvariantMass();
      pho_conv_paircotthetasep[pho_n]=conv->pairCotThetaSeparation();
      pho_conv_eoverp[pho_n]=conv->EoverP();
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
