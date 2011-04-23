#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollection.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"

#include "DataFormats/Scalers/interface/DcsStatus.h"
#include <iostream>

float GlobeElectrons::hoeCalculator(const reco::BasicCluster* clus, const CaloGeometry& geometry,
                                    const edm::Event& e , const edm::EventSetup& c) {
  
  float h = 0.;

  GlobalPoint pclu(clus->x(),clus->y(),clus->z());

  edm::Handle< HBHERecHitCollection > hbhe ;
  e.getByLabel("reducedHcalRecHits", "hbhereco",hbhe);
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


GlobeElectrons::GlobeElectrons(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  char a[100];
  sprintf(a, "ElectronColl_%s", nome);
  electronColl = iConfig.getParameter<edm::InputTag>(a);
  debug_level = iConfig.getParameter<int>("Debug_Level");
  doFastSim = iConfig.getParameter<bool>("doFastSim");
  doAodSim = iConfig.getParameter<bool>("doAodSim");
  trackColl = iConfig.getParameter<edm::InputTag>("TrackColl");
  trackColl2 = iConfig.getParameter<edm::InputTag>("TrackColl3");
  vertexColl = iConfig.getParameter<edm::InputTag>("VertexColl_std");

  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");
  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");

  eIDLabels = iConfig.getParameter<std::vector<edm::InputTag> >("eIDLabels");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeElectrons::defineBranch(TTree* tree) {

  el_sc = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_p4 = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_momvtx = new TClonesArray("TVector3", MAX_ELECTRONS);  
  el_momvtxconst = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momcalo = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momout = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_posvtx = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_poscalo = new TClonesArray("TVector3", MAX_ELECTRONS);

  el_catbased = new std::vector<std::vector<int> >;

  char a1[50], a2[50];
  
  sprintf(a1, "el_%s_n", nome);
  sprintf(a2, "el_%s_n/I", nome);
  tree->Branch(a1, &el_n, a2);
  
  sprintf(a1, "el_%s_sc", nome);
  tree->Branch(a1, "TClonesArray", &el_sc, 32000, 0);
  
  sprintf(a1, "el_%s_p4", nome);
  tree->Branch(a1, "TClonesArray", &el_p4, 32000, 0);

  sprintf(a1, "el_%s_momvtx", nome);
  tree->Branch(a1, "TClonesArray", &el_momvtx, 32000, 0);

  sprintf(a1, "el_%s_momvtxconst", nome);
  tree->Branch(a1, "TClonesArray", &el_momvtxconst, 32000, 0);

  sprintf(a1, "el_%s_momcalo", nome);
  tree->Branch(a1, "TClonesArray", &el_momcalo, 32000, 0);

  sprintf(a1, "el_%s_momout", nome);
  tree->Branch(a1, "TClonesArray", &el_momout, 32000, 0);
  
  sprintf(a1, "el_%s_posvtx", nome);
  tree->Branch(a1, "TClonesArray", &el_posvtx, 32000, 0);

  sprintf(a1, "el_%s_poscalo", nome);
  tree->Branch(a1, "TClonesArray", &el_poscalo, 32000, 0);
  
  sprintf(a1, "el_%s_eopin", nome);
  sprintf(a2, "el_%s_eopin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eopin, a2);
  
  sprintf(a1, "el_%s_eseedopout", nome);
  sprintf(a2, "el_%s_eseedopout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseedopout, a2);
  
  sprintf(a1, "el_%s_pout", nome);
  sprintf(a2, "el_%s_pout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pout, a2);
  
  sprintf(a1, "el_%s_pin", nome);
  sprintf(a2, "el_%s_pin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pin, a2);
  
  sprintf(a1, "el_%s_e1x5", nome);
  sprintf(a2, "el_%s_e1x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e1x5, a2);
  
  sprintf(a1, "el_%s_e5x5", nome);
  sprintf(a2, "el_%s_e5x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e5x5, a2);

  sprintf(a1, "el_%s_e2x5", nome);
  sprintf(a2, "el_%s_e2x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e2x5, a2);
 
  sprintf(a1, "el_%s_spp", nome);
  sprintf(a2, "el_%s_spp[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_spp, a2);
  
  sprintf(a1, "el_%s_see", nome);
  sprintf(a2, "el_%s_see[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_see, a2);

  sprintf(a1, "el_%s_sieie", nome);
  sprintf(a2, "el_%s_sieie[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_sieie, a2);

  sprintf(a1, "el_%s_sieiesc", nome);
  sprintf(a2, "el_%s_sieiesc[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_sieiesc, a2);

  sprintf(a1, "el_%s_eseed", nome);
  sprintf(a2, "el_%s_eseed[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseed, a2);

  sprintf(a1, "el_%s_enearbcopin", nome);
  sprintf(a2, "el_%s_enearbcopin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_enearbcopin, a2);

  sprintf(a1, "el_%s_eseedopin", nome);
  sprintf(a2, "el_%s_eseedopin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseedopin, a2);
  
  sprintf(a1, "el_%s_fbrem", nome);
  sprintf(a2, "el_%s_fbrem[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_fbrem, a2);

  sprintf(a1, "el_%s_nbrem", nome);
  sprintf(a2, "el_%s_nbrem[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_nbrem, a2);

  sprintf(a1, "el_%s_h", nome);
  sprintf(a2, "el_%s_h[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_h, a2);
  
  sprintf(a1, "el_%s_hoe", nome);
  sprintf(a2, "el_%s_hoe[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoe, a2);

  sprintf(a1, "el_%s_hoed1", nome);
  sprintf(a2, "el_%s_hoed1[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoed1, a2);

  sprintf(a1, "el_%s_hoed2", nome);
  sprintf(a2, "el_%s_hoed2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoed2, a2);
  
  sprintf(a1, "el_%s_detain", nome);
  sprintf(a2, "el_%s_detain[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_detain, a2);
  
  sprintf(a1, "el_%s_dphiin", nome);
  sprintf(a2, "el_%s_dphiin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dphiin, a2);
  
  sprintf(a1, "el_%s_detaout", nome);
  sprintf(a2, "el_%s_detaout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_detaout, a2);
  
  sprintf(a1, "el_%s_dphiout", nome);
  sprintf(a2, "el_%s_dphiout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dphiout, a2);
   
  sprintf(a1, "el_%s_class", nome);
  sprintf(a2, "el_%s_class[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_class, a2);
 
  sprintf(a1, "el_%s_crack", nome);
  sprintf(a2, "el_%s_crack[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_crack, a2);
   
  sprintf(a1, "el_%s_nambtk", nome);
  sprintf(a2, "el_%s_nambtk[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_nambtk, a2);

  sprintf(a1, "el_%s_scind", nome);
  sprintf(a2, "el_%s_scind[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_scind, a2);
  
  sprintf(a1, "el_%s_z0", nome);
  sprintf(a2, "el_%s_z0[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_z0, a2);
  
  sprintf(a1, "el_%s_d0", nome);
  sprintf(a2, "el_%s_d0[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_d0, a2);

  sprintf(a1, "el_%s_qoverperr", nome);
  sprintf(a2, "el_%s_qoverperr[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_qoverperr, a2);
  
  sprintf(a1, "el_%s_pterr", nome);
  sprintf(a2, "el_%s_pterr[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pterr, a2);
  
  sprintf(a1, "el_%s_etaerr", nome);
  sprintf(a2, "el_%s_etaerr[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_etaerr, a2);
  
  sprintf(a1, "el_%s_phierr", nome);
  sprintf(a2, "el_%s_phierr[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_phierr, a2);

  sprintf(a1, "el_%s_z0err", nome);
  sprintf(a2, "el_%s_z0err[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_z0err, a2);
  
  sprintf(a1, "el_%s_d0err", nome);
  sprintf(a2, "el_%s_d0err[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_d0err, a2);
  
  sprintf(a1, "el_%s_chi2", nome);
  sprintf(a2, "el_%s_chi2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_chi2, a2);
  
  sprintf(a1, "el_%s_dof", nome);
  sprintf(a2, "el_%s_dof[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dof, a2);

  sprintf(a1, "el_%s_mva", nome);
  sprintf(a2, "el_%s_mva[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_mva, a2);
  
  sprintf(a1, "el_%s_ch_gsf", nome);
  sprintf(a2, "el_%s_ch_gsf[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ch_gsf, a2);

  sprintf(a1, "el_%s_ch_scpix", nome);
  sprintf(a2, "el_%s_ch_scpix[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ch_scpix, a2);

  sprintf(a1, "el_%s_charge", nome);
  sprintf(a2, "el_%s_charge[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_charge, a2);
  
  sprintf(a1, "el_%s_losthits", nome);
  sprintf(a2, "el_%s_losthits[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_losthits, a2);
  
  sprintf(a1, "el_%s_validhits", nome);
  sprintf(a2, "el_%s_validhits[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_validhits, a2);

  sprintf(a1, "el_%s_hp_expin", nome);
  sprintf(a2, "el_%s_hp_expin[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_hp_expin, a2);

  sprintf(a1, "el_%s_hp_expout", nome);
  sprintf(a2, "el_%s_hp_expout[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_hp_expout, a2);

  sprintf(a1, "el_%s_rohighe", nome);
  sprintf(a2, "el_%s_rohighe[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_rohighe, a2);
  
  sprintf(a1, "el_%s_roloose", nome);
  sprintf(a2, "el_%s_roloose[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_roloose, a2);
  
  sprintf(a1, "el_%s_rotight", nome);
  sprintf(a2, "el_%s_rotight[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_rotight, a2);
  
  sprintf(a1, "el_%s_loose", nome);
  sprintf(a2, "el_%s_loose[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_loose, a2);
  
  sprintf(a1, "el_%s_tight", nome);
  sprintf(a2, "el_%s_tight[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_tight, a2);

  sprintf(a1, "el_%s_catbased", nome);
  tree->Branch(a1, "std::vector<std::vector<int> >", &el_catbased);

  sprintf(a1, "el_%s_tkind", nome);
  sprintf(a2, "el_%s_tkind[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_tkind, a2); 

  sprintf(a1, "el_%s_hcaliso03", nome);
  sprintf(a2, "el_%s_hcaliso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hcaliso03, a2);

  sprintf(a1, "el_%s_ecaliso03", nome);
  sprintf(a2, "el_%s_ecaliso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ecaliso03, a2);
  
  sprintf(a1, "el_%s_tkiso03", nome);
  sprintf(a2, "el_%s_tkiso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_tkiso03, a2);

  sprintf(a1, "el_%s_hcaliso04", nome);
  sprintf(a2, "el_%s_hcaliso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hcaliso04, a2);

  sprintf(a1, "el_%s_ecaliso04", nome);
  sprintf(a2, "el_%s_ecaliso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ecaliso04, a2);
  
  sprintf(a1, "el_%s_tkiso04", nome);
  sprintf(a2, "el_%s_tkiso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_tkiso04, a2);
  
  sprintf(a1, "el_%s_tkdrv", nome);
  sprintf(a2, "el_%s_tkdrv[el_%s_n]/O", nome, nome);
  tree->Branch(a1, &el_tkdrv, a2);
 
  sprintf(a1, "el_%s_ecaldrv", nome);
  sprintf(a2, "el_%s_ecaldrv[el_%s_n]/O", nome, nome);
  tree->Branch(a1, &el_ecaldrv, a2);

  sprintf(a1, "el_%s_ip_ctf", nome);
  sprintf(a2, "el_%s_ip_ctf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ip_ctf, a2);

  sprintf(a1, "el_%s_ip_gsf", nome);
  sprintf(a2, "el_%s_ip_gsf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ip_gsf, a2);

  sprintf(a1, "el_%s_dist", nome);
  sprintf(a2, "el_%s_dist[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dist, a2);

  sprintf(a1, "el_%s_dcot", nome);
  sprintf(a2, "el_%s_dcot[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dcot, a2);

  sprintf(a1, "el_%s_hp_1pxb", nome);
  sprintf(a2, "el_%s_hp_1pxb[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_1pxb, a2);

  sprintf(a1, "el_%s_hp_1pxf", nome);
  sprintf(a2, "el_%s_hp_1pxf[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_1pxf, a2);
}


bool GlobeElectrons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(electronColl, elH);
  
  edm::Handle<reco::TrackCollection> tkH;
  iEvent.getByLabel(trackColl, tkH);

  edm::Handle<reco::GsfTrackCollection> tkH2;
  iEvent.getByLabel(trackColl2, tkH2);
  
  edm::Handle<reco::SuperClusterCollection> superClustersBarrelH; 
  iEvent.getByLabel(hybridSuperClusterColl,superClustersBarrelH);
  
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vertexColl, vtxH);

  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();
  
  edm::Handle<EERecHitCollection> pEERecHitH;
  edm::Handle<EBRecHitCollection> pEBRecHitH;
  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  const EcalRecHitCollection *barrelRecHits = pEBRecHitH.product();
  const EcalRecHitCollection *endcapRecHits = pEERecHitH.product();
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;

  //Read eID results
  std::vector<edm::Handle<edm::ValueMap<float> > > eIDValueMap(5); 
  //Robust-Loose 
  iEvent.getByLabel("eidRobustLoose", eIDValueMap[0]); 
  const edm::ValueMap<float> & eIDmapRL = *eIDValueMap[0];
  //Robust-Tight 
  iEvent.getByLabel("eidRobustTight", eIDValueMap[1]); 
  const edm::ValueMap<float> & eIDmapRT = *eIDValueMap[1];
  //Loose 
  iEvent.getByLabel("eidLoose", eIDValueMap[2]); 
  const edm::ValueMap<float> & eIDmapL = *eIDValueMap[2];
  //Tight 
  iEvent.getByLabel("eidTight", eIDValueMap[3]); 
  const edm::ValueMap<float> & eIDmapT = *eIDValueMap[3];
  iEvent.getByLabel("eidRobustHighEnergy", eIDValueMap[4]); 
  const edm::ValueMap<float> & eIDmapRHE = *eIDValueMap[4];
  
  el_sc->Clear();
  el_p4->Clear(); 
  el_momvtx->Clear();
  el_momvtxconst->Clear();
  el_momcalo->Clear();
  el_posvtx->Clear();
  el_poscalo->Clear(); 
  el_catbased ->clear();
  el_n = 0;
  
  if (debug_level > 9)
    std::cout << "GlobeElectrons: Electron collection size: "<< elH->size() << std::endl;
  
  for(reco::GsfElectronCollection::const_iterator igsf = elH->begin(); igsf != elH->end(); igsf++) {

    if (el_n >= MAX_ELECTRONS) {
      std::cout << "GlobeElectrons: WARNING TOO MANY ELECTRONS: " << elH->size() << " (allowed " << MAX_ELECTRONS << ")" << std::endl;
      break;
    }
    
    reco::GsfElectron egsf = reco::GsfElectron(*igsf);
    
    if(gCUT->cut(egsf))
      continue;
    
    float phi = egsf.superCluster()->phi();
    float theta = (2*atan(exp(-egsf.superCluster()->eta())));
    float en = egsf.superCluster()->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);

    new ((*el_sc)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_sc->At(el_n))->SetXYZT(px, py, pz, en);
    
    new ((*el_p4)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_p4->At(el_n))->SetXYZT(egsf.px(), egsf.py(), egsf.pz(), egsf.energy());

    new ((*el_momvtxconst)[el_n]) TVector3();
    ((TVector3 *)el_momvtxconst->At(el_n))->SetXYZ(egsf.trackMomentumAtVtxWithConstraint().x(), 
                                              egsf.trackMomentumAtVtxWithConstraint().y(), egsf.trackMomentumAtVtxWithConstraint().z());
    
    new ((*el_momvtx)[el_n]) TVector3();
    ((TVector3 *)el_momvtx->At(el_n))->SetXYZ(egsf.trackMomentumAtVtx().x(), 
                                              egsf.trackMomentumAtVtx().y(), egsf.trackMomentumAtVtx().z());
    
    new ((*el_momcalo)[el_n]) TVector3();
    ((TVector3 *)el_momcalo->At(el_n))->SetXYZ(egsf.trackMomentumAtCalo().x(), 
                                               egsf.trackMomentumAtCalo().y(), egsf.trackMomentumAtCalo().z());

    new ((*el_momout)[el_n]) TVector3();
    ((TVector3 *)el_momout->At(el_n))->SetXYZ(egsf.trackMomentumOut().x(), 
                                               egsf.trackMomentumOut().y(), egsf.trackMomentumOut().z());
    
    new ((*el_posvtx)[el_n]) TVector3();
    ((TVector3 *)el_posvtx->At(el_n))->SetXYZ(egsf.trackPositionAtVtx().x(), 
                                              egsf.trackPositionAtVtx().y(), egsf.trackPositionAtVtx().z());

    new ((*el_poscalo)[el_n]) TVector3();
    ((TVector3 *)el_poscalo->At(el_n))->SetXYZ(egsf.trackPositionAtCalo().x(), 
                                               egsf.trackPositionAtCalo().y(), egsf.trackPositionAtCalo().z());
    
    el_pout[el_n] = egsf.trackMomentumOut().R();
    el_pin[el_n] = egsf.trackMomentumAtVtx().R();
    
    el_1pxb[el_n] = 0;
    el_1pxf[el_n] = 0;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel())
      el_1pxb[el_n] = 1;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelEndcap())
      el_1pxf[el_n] = 1;

    el_fbrem[el_n] = egsf.fbrem();

    el_h[el_n] = hoeCalculator(&(*(egsf.superCluster()->seed())), geometry, iEvent, iSetup);

    el_eseed[el_n] = egsf.superCluster()->seed()->energy();
    el_eseedopout[el_n] = egsf.eSeedClusterOverPout();
    el_eseedopin[el_n] = egsf.eSeedClusterOverP();
    el_enearbcopin[el_n] = egsf.electronCluster()->energy()/egsf.trackMomentumAtVtx().R();
    el_eopin[el_n] = egsf.eSuperClusterOverP();
    el_detain[el_n] = egsf.deltaEtaSuperClusterTrackAtVtx();
    el_dphiin[el_n] = egsf.deltaPhiSuperClusterTrackAtVtx();
    el_detaout[el_n] = egsf.deltaEtaSeedClusterTrackAtCalo();
    el_dphiout[el_n] = egsf.deltaPhiSeedClusterTrackAtCalo();

    el_nambtk[el_n] = egsf.ambiguousGsfTracksSize();
    el_class[el_n] = egsf.classification();

    el_nbrem[el_n] = egsf.numberOfBrems();
   
    el_e5x5[el_n] = egsf.e5x5();
    el_e2x5[el_n] = egsf.e2x5Max();
    el_e1x5[el_n] = egsf.e1x5();

    el_sieie[el_n] = egsf.sigmaIetaIeta();
    el_see[el_n] = egsf.sigmaEtaEta();

    el_hoe[el_n] = egsf.hcalOverEcal();
    el_hoed1[el_n] = egsf.hcalDepth1OverEcal();
    el_hoed2[el_n] = egsf.hcalDepth2OverEcal();

    el_d0[el_n] = egsf.gsfTrack()->d0();
    el_z0[el_n] = egsf.gsfTrack()->dz(); 
    el_qoverperr[el_n] = egsf.gsfTrack()->qoverpError();
    el_pterr[el_n] = egsf.gsfTrack()->ptError();
    el_etaerr[el_n] = egsf.gsfTrack()->etaError();
    el_phierr[el_n] = egsf.gsfTrack()->phiError();
    el_d0err[el_n] = egsf.gsfTrack()->d0Error();
    el_z0err[el_n] = egsf.gsfTrack()->dzError();
    el_chi2[el_n] = egsf.gsfTrack()->chi2();
    el_dof[el_n] = egsf.gsfTrack()->ndof();
    el_validhits[el_n] = egsf.gsfTrack()->numberOfValidHits();
    el_losthits[el_n] = egsf.gsfTrack()->numberOfLostHits();

    math::XYZPoint vtxPoint(0.0,0.0,0.0);
    if (vtxH->size() != 0) {
      reco::VertexRef vtx(vtxH, 0);
      vtxPoint = math::XYZPoint(vtx->x(),vtx->y(),vtx->z());
    }

    el_ip_gsf[el_n] = egsf.gsfTrack()->dxy(vtxPoint);

    if (!doAodSim && (trackColl2.encode() != "electronGsfTracks")) {
      for(unsigned int j=0; j<tkH2->size(); j++) { 
        reco::GsfTrackRef tk2(tkH2, j);
        std::pair<unsigned int, float> result = sharedHits(*tk2, *(egsf.gsfTrack()));
        
        if (result.second > .999) {
          el_hp_expin[el_n] = tk2->trackerExpectedHitsInner().numberOfHits();
          el_hp_expout[el_n] = tk2->trackerExpectedHitsOuter().numberOfHits();
          break;
        }
      }
    } else {
      el_hp_expin[el_n] = egsf.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      el_hp_expout[el_n] = egsf.gsfTrack()->trackerExpectedHitsOuter().numberOfHits();
    }
 
    el_charge[el_n] = egsf.charge();
    
    el_ch_scpix[el_n] = egsf.scPixCharge();
    el_ch_gsf[el_n] = egsf.gsfTrack()->charge();

    if (egsf.isEB())
      el_crack[el_n] = 0;

    if (egsf.isEE())
      el_crack[el_n] = 1;

    if (egsf.isEBEEGap())
      el_crack[el_n] = 3;

    if (egsf.isEBGap())
      el_crack[el_n] = 4;

    if (egsf.isEEGap())
      el_crack[el_n] = 5;

    int cmsTkind = -1;

    el_tkind[el_n] = -1;
    el_scind[el_n] = -1;
  
    if(egsf.closestCtfTrackRef().isNonnull()) {
      for(unsigned int j=0; j<tkH->size(); ++j) {
        reco::TrackRef tk(tkH, j);
        if(gCUT->cut(*tk))
          continue; 
        if (tk == egsf.closestCtfTrackRef()) {
          el_tkind[el_n] = j;
          cmsTkind = j;
          el_ip_ctf[el_n] = egsf.closestCtfTrackRef()->dxy(vtxPoint);
        }
      }
    } else {
      el_tkind[el_n] = -1;
      el_ip_ctf[el_n] = 0;
    }
    
    int index = 0;
    // loop over the two SC collections
    for(int z = 0; z<2; ++z) {
      if (z == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersBarrelH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersBarrelH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))continue;
          // passed cuts
          
          if (&(*egsf.superCluster()) == &(*cluster)) {
            el_scind[el_n] = index; 
            el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*barrelRecHits), &(*topology))[0]);
            std::vector<float> vCov = EcalClusterTools::covariances( *(cluster->seed()), &(*barrelRecHits), &(*topology), &geometry);
            el_spp[el_n] = sqrt(vCov[2]);
            break;
          }
          index++;
        }
      }
      
      if (z == 1) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersEndcapH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))continue;
          // passed cuts
          
          if (&(*(egsf.superCluster())) == &(*cluster)) {
            el_scind[el_n] = index;
            el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*endcapRecHits), &(*topology))[0]);
            std::vector<float> vCov = EcalClusterTools::covariances( *(cluster->seed()), &(*endcapRecHits), &(*topology), &geometry);
            el_spp[el_n] = sqrt(vCov[2]);
            break;
          }
          index++;
        }
      }
    }

    if (el_scind[el_n] == -1) {
      if(fabs(egsf.superCluster()->eta()) < 1.479) {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*barrelRecHits), &(*topology))[0]);
      } else {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*endcapRecHits), &(*topology))[0]);
      }
    }
    
    //electronID
    reco::GsfElectronRef electronRef(elH, std::distance(elH->begin(), igsf));
    el_roloose[el_n] = (Int_t)eIDmapRL[electronRef];
    el_rotight[el_n] = (Int_t)eIDmapRT[electronRef];
    el_rohighe[el_n] = (Int_t)eIDmapRHE[electronRef];
    el_loose[el_n] = (Int_t)eIDmapL[electronRef];
    el_tight[el_n] = (Int_t)eIDmapT[electronRef];
    
    el_mva[el_n] = egsf.mva();

    el_ecaldrv[el_n] = egsf.ecalDrivenSeed();
    el_tkdrv[el_n] = egsf.trackerDrivenSeed();
    
    el_tkiso04[el_n] = egsf.dr04TkSumPt();
    el_ecaliso04[el_n] = egsf.dr04EcalRecHitSumEt();
    el_hcaliso04[el_n] = egsf.dr04HcalTowerSumEt();

    el_tkiso03[el_n] = egsf.dr03TkSumPt();
    el_ecaliso03[el_n] = egsf.dr03EcalRecHitSumEt();
    el_hcaliso03[el_n] = egsf.dr03HcalTowerSumEt();

    // Fill out electron identification
    std::vector<edm::Handle<edm::ValueMap<float> > > eIDVM(9); 
    std::vector<int> results;
    
    for(unsigned int j=0; j<eIDLabels.size(); j++) {
      if (iEvent.getByLabel(eIDLabels[j], eIDVM[j])) {  
        const edm::ValueMap<float>& eIDmapTemp = *eIDVM[j];
        reco::GsfElectronRef electronRef(elH, std::distance(elH->begin(), igsf));
        results.push_back((Int_t)eIDmapTemp[electronRef]);
      } else {
        results.push_back(-1);
      }
    }
    
    el_catbased->push_back(results);

    el_dist[el_n] = (egsf.convDist() == -9999.? 9999:egsf.convDist());
    el_dcot[el_n] = (egsf.convDcot() == -9999.? 9999:egsf.convDcot());
    
    el_n++;
  }
  
  return true;
}


bool GlobeElectrons::inCrack(float etain) {

  float eta = fabs(etain);
  
  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

std::pair<unsigned int, float> GlobeElectrons::sharedHits(const reco::Track& trackA, const reco::Track& trackB) {
  
  unsigned int shared = 0;
  for(trackingRecHit_iterator tkHitA = trackA.recHitsBegin(); tkHitA !=trackA.recHitsEnd(); ++tkHitA){
    for(trackingRecHit_iterator tkHitB = trackB.recHitsBegin();
        tkHitB !=trackB.recHitsEnd(); ++tkHitB){
      if( (**tkHitA).isValid() && (**tkHitB).isValid() &&(**tkHitA).sharesInput( &(**tkHitB),TrackingRecHit::all)) {
        shared++;
        break;
      }
    }
  }

  float fraction = (float) shared/std::min(trackA.found(),trackB.found());
  return std::make_pair(shared,fraction);
}

/*
bool GlobeElectrons::analyze_pf(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Candidate info
  edm::Handle<reco::PFCandidateCollection> collection;
  edm::InputTag label("particleFlow","");
  iEvent.getByLabel(label, collection);
  std::vector<reco::PFCandidate> candidates = (*collection.product());
  std::vector<reco::PFCandidate>::iterator it;
  
  el_p4->Clear(); 
  el_momvtx->Clear();
  el_momcalo->Clear();
  el_posvtx->Clear();
  el_n = 0;

  for (it = candidates.begin(); it != candidates.end(); ++it )   {
    reco::PFCandidate::ParticleType type = (*it).particleId();
    if (type == reco::PFCandidate::e)    {
      // in PFCandidates the chosen cut is mva > -0.1

      if (el_n >= MAX_ELECTRONS) {
        std::cout << "GlobeElectrons: WARNING TOO MANY ELECTRONS: " << candidates.size() << " (allowed " << MAX_ELECTRONS << ")" << std::endl;
        break;
      }
    
      reco::PFCandidate egsf = reco::PFCandidate(*it);
      new ((*el_p4)[el_n]) TLorentzVector();
      ((TLorentzVector *)el_p4->At(el_n))->SetXYZT(egsf.px(), egsf.py(), egsf.pz(), egsf.energy());
      
      new ((*el_momvtx)[el_n]) TVector3();
      ((TVector3 *)el_momvtx->At(el_n))->SetXYZ(0,0,0);
      
      new ((*el_momcalo)[el_n]) TVector3();
      ((TVector3 *)el_momcalo->At(el_n))->SetXYZ(0,0,0);      
      new ((*el_posvtx)[el_n]) TVector3();
      ((TVector3 *)el_momcalo->At(el_n))->SetXYZ(0,0,0);
      
      el_pout[el_n] = 0;
      el_pin[el_n] = 0;
      el_eseedopout[el_n] = 0;
      el_eopin[el_n] = 0;
      el_fbrem[el_n] = 0;
      
      el_sieie[el_n] = 0;
      el_see[el_n] = 0;
      
      el_eseed[el_n] = 0;
      el_eseedopin[el_n] = 0;
      el_hoe[el_n] = 0;
      el_detain[el_n] = 0;
      el_dphiin[el_n] = 0;
      el_detaout[el_n] = 0;
      el_dphiout[el_n] = 0;
      el_class[el_n] = 0;
      el_d0[el_n] = 0;
      el_z0[el_n] = 0;
      el_qoverperr[el_n] = 0;
      el_pterr[el_n] = 0;
      el_etaerr[el_n] = 0;
      el_phierr[el_n] = 0;
      el_d0err[el_n] = 0;
      el_z0err[el_n] = 0;
      el_chi2[el_n] = 0;
      el_dof[el_n] = 0;
      el_validhits[el_n] = 0;
      el_losthits[el_n] = 0;
      el_charge[el_n] = 0;
      el_crack[el_n] = 0;
      

      el_tkind[el_n] = -1;
      el_scind[el_n] = -1;
      
      el_roloose[el_n] = -1;
      el_rotight[el_n] = -1;
      el_rohighe[el_n] = -1;
      el_loose[el_n] = -1;
      el_tight[el_n] = -1;
      
      el_mva[el_n] = egsf.mva_e_pi();
      
      el_tkdrv[el_n] = 0;
      el_ecaldrv[el_n] = 0;
    
      el_tkiso03[el_n] = 0;
      el_ecaliso03[el_n] = 0;
      el_hcaliso03[el_n] = 0;

      el_tkiso04[el_n] = 0;
      el_ecaliso04[el_n] = 0;
      el_hcaliso04[el_n] = 0;


      //el_3dip_valid[el_n] = false;
      //el_3dip_x[el_n]     = -9999;
      //el_3dip_y[el_n]     = -9999;
      //el_3dip_z[el_n]     = -9999;
      //el_3dip_xerr[el_n]  = -9999;
      //el_3dip_yerr[el_n]  = -9999;
      //el_3dip_zerr[el_n]  = -9999;
      el_n++;
    }
  }
  
  return true;
}
*/

