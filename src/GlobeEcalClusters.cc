#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include <vector>
#include <algorithm>

//----------------------------------------------------------------------

GlobeEcalClusters::GlobeEcalClusters(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  debug_level = iConfig.getParameter<int>("Debug_Level");
  gCUT = new GlobeCuts(iConfig);
  
  edm::ParameterSet psetSC = iConfig.getParameter<edm::ParameterSet>("SuperClusterCuts");
  edm::ParameterSet psetBC = iConfig.getParameter<edm::ParameterSet>("BasicClusterCuts");
  
  //super clusters
  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");
  
  //basic clusters
  barrelHybridClusterColl = iConfig.getParameter<edm::InputTag>("BarrelHybridClusterColl");
  barrelBasicClusterColl = iConfig.getParameter<edm::InputTag>("BarrelBasicClusterColl");
  endcapBasicClusterColl = iConfig.getParameter<edm::InputTag>("EndcapBasicClusterColl");
  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");
}

//----------------------------------------------------------------------

void GlobeEcalClusters::defineBranch(TTree* tree) {
  char a2[100];
  sc_p4 = new TClonesArray("TLorentzVector", MAX_SUPERCLUSTERS);
  sc_xyz = new TClonesArray("TVector3", MAX_SUPERCLUSTERS);
  sc_islbar_p4 = new TClonesArray("TLorentzVector", MAX_SUPERCLUSTERS);
  sc_islbar_xyz = new TClonesArray("TVector3", MAX_SUPERCLUSTERS);
  bc_p4 = new TClonesArray("TLorentzVector", MAX_BASICCLUSTERS);
  bc_xyz = new TClonesArray("TVector3", MAX_SUPERCLUSTERS);
  
  tree->Branch("sc_islbar_n", &sc_islbar_n, "sc_islbar_n/I");
  tree->Branch("sc_islbar_p4", "TClonesArray", &sc_islbar_p4, 32000, 0);
  tree->Branch("sc_islbar_xyz", "TClonesArray", &sc_islbar_xyz, 32000, 0);
  tree->Branch("sc_islbar_raw", &sc_islbar_raw, "sc_islbar_raw[sc_islbar_n]/F");
  tree->Branch("sc_islbar_nbc", &sc_islbar_nbc, "sc_islbar_nbc[sc_islbar_n]/I");
  tree->Branch("sc_islbar_seedenergy", &sc_islbar_seedenergy,"sc_islbar_seedenergy[sc_islbar_n]/F" );
  tree->Branch("sc_islbar_bcseedind", &sc_islbar_bcseedind, "sc_islbar_bcseedind[sc_islbar_n]/I");
  sprintf (a2, "sc_islbar_bcind[sc_islbar_n][%d]/I", MAX_SUPERCLUSTER_BASICCLUSTERS);
  tree->Branch("sc_islbar_bcind", &sc_islbar_bcind, a2);
  
  //SC hybrid in barrel and island in endcap
  tree->Branch("sc_n", &sc_n, "sc_n/I");
  tree->Branch("sc_islend_n", &sc_islend_n, "sc_islend_n/I");
  tree->Branch("sc_hybrid_n", &sc_hybrid_n, "sc_hybrid_n/I");
  tree->Branch("sc_p4", "TClonesArray", &sc_p4, 32000, 0);
  tree->Branch("sc_xyz", "TClonesArray", &sc_xyz, 32000, 0);
  tree->Branch("sc_pre", &sc_pre, "sc_pre[sc_n]/F");
  tree->Branch("sc_raw", &sc_raw, "sc_raw[sc_n]/F");
  tree->Branch("sc_barrel", &sc_barrel, "sc_barrel[sc_n]/I");
  tree->Branch("sc_2xN", &sc_2xN, "sc_2xN[sc_n]/F");
  tree->Branch("sc_5xN", &sc_5xN, "sc_5xN[sc_n]/F");
  tree->Branch("sc_sieie", &sc_sieie, "sc_sieie[sc_n]/F");
#ifdef CMSSW_VERSION_210
  tree->Branch("sc_see", &sc_see, "sc_see[sc_n]/F");
#endif
  tree->Branch("sc_nbc", &sc_nbc, "sc_nbc[sc_n]/I");
  tree->Branch("sc_bcseedind", &sc_bcseedind, "sc_bcseedind[sc_n]/I");
  sprintf (a2, "sc_bcind[sc_n][%d]/I", MAX_SUPERCLUSTER_BASICCLUSTERS);
  tree->Branch("sc_bcind", &sc_bcind, a2);
  
  //CHECK add shape variables for hybrid
  //CHECK  InputTag BarrelHybridClusterShapeColl = hybridSuperClusters:hybridShapeAssoc
  
  //basic clusters
  tree->Branch("bc_n", &bc_n, "bc_n/I");
  tree->Branch("bc_islbar_n", &bc_islbar_n, "bc_islbar_n/I");
  tree->Branch("bc_islend_n", &bc_islend_n, "bc_islend_n/I");
  tree->Branch("bc_hybrid_n", &bc_hybrid_n, "bc_hybrid_n/I");
  tree->Branch("bc_p4", "TClonesArray", &bc_p4, 32000, 0);
  tree->Branch("bc_xyz", "TClonesArray", &bc_xyz, 32000, 0);
  tree->Branch("bc_nhits", &bc_nhits,"bc_nhits[bc_n]/I");
  tree->Branch("bc_s1", &bc_s1, "bc_s1[bc_n]/F");
  tree->Branch("bc_rook", &bc_rook, "bc_rook[bc_n]/F");
  tree->Branch("bc_s4", &bc_s4, "bc_s4[bc_n]/F");
  tree->Branch("bc_s9", &bc_s9, "bc_s9[bc_n]/F");
  tree->Branch("bc_s25", &bc_s25, "bc_s25[bc_n]/F");
  //tree->Branch("bc_hoe", &bc_hoe, "bc_hoe[bc_n]/F");
  //tree->Branch("bc_radius", &bc_radius, "bc_radius[bc_n]/F");
  //tree->Branch("bc_z", &bc_z, "bc_z[bc_n]/F");
  tree->Branch("bc_spp", &bc_spp, "bc_spp[bc_n]/F");
  tree->Branch("bc_see", &bc_see, "bc_see[bc_n]/F");
  tree->Branch("bc_sep", &bc_sep, "bc_sep[bc_n]/F");
  tree->Branch("bc_type", &bc_type, "bc_type[bc_n]/I");//type 1 = hybrid, 2 = island endcap, 3 = island barrel.

  tree->Branch("bc_s1x5_0", &bc_s1x5_0, "bc_s1x5_0[bc_n]/F");
  tree->Branch("bc_s1x5_1", &bc_s1x5_1, "bc_s1x5_1[bc_n]/F");
  tree->Branch("bc_s1x5_2", &bc_s1x5_2, "bc_s1x5_2[bc_n]/F");
  tree->Branch("bc_s1x3_0", &bc_s1x3_0, "bc_s1x3_0[bc_n]/F");
  tree->Branch("bc_s1x3_1", &bc_s1x3_1, "bc_s1x3_1[bc_n]/F");
  tree->Branch("bc_s1x3_2", &bc_s1x3_2, "bc_s1x3_2[bc_n]/F");
  tree->Branch("bc_s3x1_0", &bc_s3x1_0, "bc_s3x1_0[bc_n]/F");
  tree->Branch("bc_s3x1_1", &bc_s3x1_1, "bc_s3x1_1[bc_n]/F");
  tree->Branch("bc_s3x1_2", &bc_s3x1_2, "bc_s3x1_2[bc_n]/F");
  tree->Branch("bc_s5x1_0", &bc_s5x1_0, "bc_s5x1_0[bc_n]/F");
  tree->Branch("bc_s5x1_1", &bc_s5x1_1, "bc_s5x1_1[bc_n]/F");
  tree->Branch("bc_s5x1_2", &bc_s5x1_2, "bc_s5x1_2[bc_n]/F");

  tree->Branch("bc_sieie", &bc_sieie, "bc_sieie[bc_n]/F");
}

//----------------------------------------------------------------------

bool GlobeEcalClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get the collection geometry:
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  geometry = geoHandle.product();

  // Get EcalRecHits
  edm::Handle<EBRecHitCollection> pEBRecHitH;
  edm::Handle<EERecHitCollection> pEERecHitH;
  //edm::Handle<ESRecHitCollection> pESRecHitH; 
  
  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  
  barrelRecHits = pEBRecHitH.product();
  endcapRecHits = pEERecHitH.product();

  edm::ESHandle<CaloTopology> theCaloTopo;
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  topology = theCaloTopo.product();
 
  // get collections

  iEvent.getByLabel(hybridSuperClusterColl,superClustersHybridH);
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  const std::string instan="hybridBarrelBasicClusters";
  iEvent.getByLabel(barrelHybridClusterColl, hybridClustersBarrelH);
  iEvent.getByLabel(endcapBasicClusterColl, basicClustersEndcapH);
 
  if (debug_level > 9) {
    std::cout << "GlobeEcalClusters: superClustersEndcapH->size() "<< superClustersEndcapH->size() << std::endl;
    std::cout << "GlobeEcalClusters: hybridClustersBarrelH->size() "<< hybridClustersBarrelH->size() << std::endl;
    std::cout << "GlobeEcalClusters: basicClustersEndcapH->size() "<< basicClustersEndcapH->size() << std::endl;
  }
  


  //----------------------------------------    
  // analyze super clusters
  //----------------------------------------  
  {
    sc_p4->Clear();
    sc_xyz->Clear();
    sc_islbar_p4->Clear();
    sc_islbar_xyz->Clear();

    sc_n = 0;
    sc_hybrid_n = 0;
    sc_islend_n = 0;
    sc_islbar_n = 0;

    // BARREL 
    analyzeBarrelSuperClusters();
    
    // ISLAND ENDCAP
    analyzeEndcapSuperClusters();
  }
  
  //----------------------------------------  
  // analyze basic clusters 
  //----------------------------------------  
  {
    bc_n = 0;
    bc_islbar_n = 0;
    bc_islend_n = 0;
    bc_hybrid_n = 0;

    bc_p4->Clear();
    bc_xyz->Clear();

    // endcap (note the order that endcap comes first)
    analyzeEndcapBasicClusters();

    // barrel
    analyzeBarrelHybridClusters();
  }

  
  //if (bc_n == 0 || sc_islbar_n == 0)
  //  return false;
  
  return true;
}

//----------------------------------------------------------------------

void GlobeEcalClusters::e2xNOe5xN(float& e2xN, float& e5xN, const reco::SuperCluster* superCluster, const EcalRecHitCollection* hits,
                                  const  CaloSubdetectorTopology* topology) {
  
  // Take hits DetId of the cluster
#ifdef CMSSW_VERSION_210
  std::vector<DetId> crystals = superCluster->getHitsByDetId();
#else
  std::vector<std::pair<DetId, float> > crystals = superCluster->hitsAndFractions();
#endif

  // look for the max energy crystal
  double eMax=0;
  DetId eMaxId(0);
#ifdef CMSSW_VERSION_210
  std::vector<DetId>::iterator it; 
#else
  std::vector<std::pair<DetId, float> >::iterator it;
#endif
  EcalRecHitCollection::const_iterator itt;
  
  EcalRecHit testEcalRecHit;
  
  for(it = crystals.begin(); it != crystals.end(); it++) {
    
#ifdef CMSSW_VERSION_210
    itt = hits->find(*it);
    if ((*it != DetId(0)) && (itt != hits->end())) {
#else
    itt = hits->find((*it).first);
    if (((*it).first != DetId(0)) && (itt != hits->end())) {
#endif
      if(itt->energy() > eMax) {
        eMax = itt->energy();
        eMaxId = itt->id();
      }
    }
  }
  
  float e2xNdown=0, e2xNup=0;
  e2xN = 0; e5xN = 0;

  CaloNavigator<DetId> posCurrent = CaloNavigator<DetId>(eMaxId, topology);

  for(unsigned int eta=0; eta<5; eta++) {
    for(unsigned int phi=0; phi<crystals.size()*2; phi++) {

      posCurrent.home();
      posCurrent.offsetBy(eta-2, phi - crystals.size()); 

      EcalRecHitCollection::const_iterator itCrystal = hits->find(*posCurrent);
      if ((*posCurrent != DetId(0)) && (itCrystal != hits->end())) {
        e5xN += itCrystal->energy();
        if (eta == 2) {
          e2xNup += itCrystal->energy();
          e2xNdown += itCrystal->energy();
        }
        
        if (eta == 1) 
          e2xNdown += itCrystal->energy();
        if (eta == 3) 
          e2xNup += itCrystal->energy();
      }
    }
  }
  
  if (e2xNdown > e2xNup) 
    e2xN = e2xNdown;
  else
    e2xN = e2xNup;

  return;
}

//----------------------------------------------------------------------
void 
GlobeEcalClusters::analyzeBarrelSuperClusters()
{
      
  for(unsigned int i=0; i<superClustersHybridH->size(); i++) {
        
    if(sc_n >= MAX_SUPERCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many superclusters. " << MAX_SUPERCLUSTERS << " allowed.  Event has " << superClustersHybridH->size() +  superClustersEndcapH->size()<< std::endl;
      break;
    }
    
    reco::SuperClusterRef sc(superClustersHybridH, i);
    // apply the cuts
    if(gCUT->cut(*sc))continue;
    // passed cuts
    
    //Old Style:
    float phi = sc->phi();
    float theta = (2*atan(exp(-sc->eta())));
    float en = sc->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);
        
    new ((*sc_p4)[sc_n]) TLorentzVector();
    ((TLorentzVector *)sc_p4->At(sc_n))->SetXYZT(px, py, pz, en);

    new ((*sc_xyz)[sc_n]) TVector3();
    ((TVector3 *)sc_xyz->At(sc_n))->SetXYZ(sc->position().x(), sc->position().y(), sc->position().z());

    sc_raw[sc_n] = sc->rawEnergy();
    //better to zero it 	 
    sc_pre[sc_n]=0;
    //e2xNOe5xN(sc_2xN[sc_n], sc_5xN[sc_n], &(*sc), barrelRecHits, topology_matteo);
    
#ifdef CMSSW_VERSION_210
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*barrelRecHits), &(*topology), &geometry)[0]);
    sc_see[sc_n] = sqrt(EcalClusterTools::teoCovariances(*(sc), &(*barrelRecHits), &(*topology), &geometry)[0]);
#else
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*barrelRecHits), &(*topology))[0]);
#endif
    
    //SEED BC 
    if (debug_level > 10)
      std::cout << sc->clustersSize() << std::endl;
    sc_nbc[sc_n]= sc->clustersSize();
        
    // get index to seed basic cluster
    for(unsigned int j=0; j<hybridClustersBarrelH->size(); ++j) {
      reco::BasicClusterRef basic(hybridClustersBarrelH, j);
      if (&(*sc->seed()) == &(*basic)) {
        sc_bcseedind[sc_n] = j + basicClustersEndcapH->size();
        break;
      }
    }
        
    // get indices to basic clusters
    if (sc->clustersSize() > 0) { 
      int limit = 0;
      
     
#ifdef CMSSW_VERSION_210
      for(reco::basicCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
#else
 for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
#endif

        if (limit >= MAX_SUPERCLUSTER_BASICCLUSTERS) {
          std::cout << "GlobeEcalCluster: WARNING too many basiclusters. in basicClustersBarrelH (" << 
            MAX_SUPERCLUSTER_BASICCLUSTERS << " allowed). Event has " << sc->clustersSize() << std::endl;
          break;
        } 
        
        for(unsigned int j=0; j<hybridClustersBarrelH->size(); ++j) {
          reco::BasicClusterRef basic(hybridClustersBarrelH, j);
          if (&(**itClus) == &(*basic)) {
            sc_bcind[sc_n][limit] = j + basicClustersEndcapH->size();
            break;
          }
        }
        limit++;
      }
    }
    sc_barrel[sc_n] = 1;
    sc_n++;
    sc_hybrid_n++;
  }
}

//----------------------------------------------------------------------  
void 
GlobeEcalClusters::analyzeEndcapSuperClusters()
{
  for(unsigned int i=0; i<superClustersEndcapH->size(); i++) {
    
    if(sc_n >= MAX_SUPERCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many superclusters. " << MAX_SUPERCLUSTERS << " allowed.  Event has " << superClustersHybridH->size() + superClustersEndcapH->size() << std::endl;
      break;
    }
    
    reco::SuperClusterRef sc(superClustersEndcapH, i);
    // apply the cuts
    if(gCUT->cut(*sc))continue;
    // passed cuts
    
    float phi = sc->phi();
    float theta = (2*atan(exp(-sc->eta())));
    float en = sc->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);
    
    new ((*sc_p4)[sc_n]) TLorentzVector();
    ((TLorentzVector *)sc_p4->At(sc_n))->SetXYZT(px, py, pz, en);
    
    new ((*sc_xyz)[sc_n]) TVector3();
    ((TVector3 *)sc_xyz->At(sc_n))->SetXYZ(sc->position().x(), sc->position().y(), sc->position().z());
    
    sc_raw[sc_n]= sc->rawEnergy();
    sc_pre[sc_n]= sc->preshowerEnergy();
    sc_nbc[sc_n]= sc->clustersSize();
    sc_2xN[sc_n] = -1; 
    sc_5xN[sc_n] = -1;
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*endcapRecHits), &(*topology))[0]);

#ifdef CMSSW_VERSION_210
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*endcapRecHits), &(*topology), &geometry)[0]);
    sc_see[sc_n] = sqrt(EcalClusterTools::teoCovariances(*(sc), &(*endcapRecHits), &(*topology), &geometry)[0]);
#else
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*endcapRecHits), &(*topology))[0]);
#endif

    // get index to seed basic cluster
    for(unsigned int j=0; j<basicClustersEndcapH->size(); ++j) {
      reco::BasicClusterRef basic(basicClustersEndcapH, j);
      if (&(*sc->seed()) == &(*basic)) {
        sc_bcseedind[sc_n] = j;
        break;
      }
    }
    
    // get indices to basic clusters
    if (sc->clustersSize() > 0) { 
      int limit = 0;
#ifdef CMSSW_VERSION_210
      for(reco::basicCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
#else
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
#endif

        if (limit >= MAX_SUPERCLUSTER_BASICCLUSTERS) {
          std::cout << "GlobeEcalCluster: WARNING too many basiclusters. in basicClustersEndcapH (" << MAX_SUPERCLUSTER_BASICCLUSTERS << " allowed). Event has " << sc->clustersSize() << std::endl;
          break;
        } 
        for(unsigned int j=0; j<basicClustersEndcapH->size(); ++j) {
          reco::BasicClusterRef basic(basicClustersEndcapH, j);
          if (&(**itClus) == &(*basic)) {
            sc_bcind[sc_n][limit] = j;
            break;
          }
        }
        limit++;
      }
    }
    
    sc_barrel[sc_n] = 0;
    sc_n++;
    sc_islend_n++;
  }
  
}


//----------------------------------------------------------------------
void 
GlobeEcalClusters::analyzeEndcapBasicClusters()
{
  // BASIC CLUSTERS ENDCAP //CHECK: missing type of clusters //MARCO
  for(unsigned int i=0; i< basicClustersEndcapH->size(); i++) {
    
    if(bc_n >= MAX_BASICCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many basicclusters. " << MAX_BASICCLUSTERS << " allowed.  Event has " << basicClustersEndcapH->size() + hybridClustersBarrelH->size() << std::endl;
      break;
    }
    
    reco::BasicClusterRef bc(basicClustersEndcapH, i);

    // make the cuts
    if(gCUT->cut(*bc))continue;
    // passed cuts
    
    float phi = bc->phi();
    float theta = (2*atan(exp(-bc->eta())));
    float en = bc->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);
    
    new ((*bc_p4)[bc_n]) TLorentzVector();
    ((TLorentzVector *)bc_p4->At(bc_n))->SetXYZT(px, py, pz, en);

    new ((*bc_xyz)[sc_n]) TVector3();
    ((TVector3 *)bc_xyz->At(sc_n))->SetXYZ(bc->position().x(), bc->position().y(), bc->position().z());

    std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
    bc_nhits[bc_n] = hits.size(); // CHECK no direct method in the dataFormat only getRecHitsByDetId

    // compute position of ECAL shower
    //float e3x3=   EcalClusterTools::e3x3(  *(bc), &(*endcapRecHits), &(*topology)); 
    //float r9 =e3x3/(aClus->rawEnergy()+aClus->preshowerEnergy());
    //float e5x5= EcalClusterTools::e5x5( *(aClus->seed()), &(*hits), &(*topology)); 
    std::pair<DetId, float> mypair=EcalClusterTools::getMaximum( *(bc), &(*endcapRecHits)); 
    
    bc_s1[bc_n] = EcalClusterTools::eMax(*(bc), &(*endcapRecHits)); 
    bc_s4[bc_n] = EcalClusterTools::e2x2(*(bc), &(*endcapRecHits), &(*topology)); 
    bc_s9[bc_n] = EcalClusterTools::e3x3(*(bc), &(*endcapRecHits), &(*topology)); 
    bc_s25[bc_n] = EcalClusterTools::e5x5(*(bc), &(*endcapRecHits), &(*topology)); 
    
    std::vector<float> rook_vect;
    rook_vect.push_back(EcalClusterTools::eLeft(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eTop(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eBottom(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eRight(*(bc), &(*endcapRecHits), &(*topology)));
    bc_rook[bc_n] = *(max_element(rook_vect.begin(), rook_vect.end()));

    std::vector<float> vCov = EcalClusterTools::covariances( *(bc), &(*endcapRecHits), &(*topology), geometry);
    
    bc_see[bc_n] = vCov[0];
    bc_sep[bc_n] = vCov[1];
    bc_spp[bc_n] = vCov[2];
    
    bc_s1x5_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), 0, 0, -2, 2);
    bc_s1x5_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -1, -1, -2, 2);
    bc_s1x5_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), 1, 1, -2, 2);
    bc_s1x3_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), 0, 0, -1, 1);
    bc_s1x3_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -1, -1, -1, 1);
    bc_s1x3_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), 1, 1, -1, 1);
    bc_s5x1_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -2, 2, 0, 0);
    bc_s5x1_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -2, 2, -1, -1);
    bc_s5x1_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -2, 2, 1, 1);
    bc_s3x1_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -1, 1, 0, 0);
    bc_s3x1_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -1, 1, -1, -1);
    bc_s3x1_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*endcapRecHits), &(*topology), mypair.first(), -1, 1, 1, 1);

    // localCovariances no longer takes geometry as input.  At this point still need to do "cvs co RecoEcal/EgammaClusterTools" from the cern cvs. REMOVE COMENT ONCE TI WORKS
    bc_sieie[bc_n] = sqrt(EcalClusterTools::localCovariances(*(bc), &(*endcapRecHits), &(*topology))[0]);
    bc_type[bc_n] = 2;

    bc_n++;
    bc_islend_n++;
  }
}

//----------------------------------------------------------------------
void
GlobeEcalClusters::analyzeBarrelHybridClusters()
{
  // HYBRID CLUSTERS BARREL //CHECK: missing type of clusters //MARCO
  for(unsigned int i=0; i< hybridClustersBarrelH->size(); i++) {
    
    if(bc_n >= MAX_BASICCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many basicclusters. " << MAX_BASICCLUSTERS << " allowed.  Event has " << basicClustersEndcapH->size() + hybridClustersBarrelH->size() << std::endl;
      break;
    }
    
    reco::BasicClusterRef bc(hybridClustersBarrelH, i);
    // make the cuts
    if(gCUT->cut(*bc))continue;
    // passed cuts
    
    float phi = bc->phi();
    float theta = (2*atan(exp(-bc->eta())));
    float en = bc->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);

    new ((*bc_p4)[bc_n]) TLorentzVector();
    ((TLorentzVector *)bc_p4->At(bc_n))->SetXYZT(px, py, pz, en);

    std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
    bc_nhits[bc_n] = hits.size(); // CHECK no direct method in the dataFormat only getRecHitsByDetId
    
    std::pair<DetId, float> mypair=EcalClusterTools::getMaximum( *(bc), &(*barrelRecHits)); 
    bc_s1[bc_n] = EcalClusterTools::eMax(  *(bc), &(*barrelRecHits));
    bc_s4[bc_n] = EcalClusterTools::e2x2(  *(bc), &(*barrelRecHits), &(*topology));
    bc_s9[bc_n] = EcalClusterTools::e3x3(  *(bc), &(*barrelRecHits), &(*topology));
    bc_s25[bc_n] = EcalClusterTools::e5x5(  *(bc), &(*barrelRecHits), &(*topology));

    std::vector<float> rook_vect;
    rook_vect.push_back(EcalClusterTools::eLeft(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eTop(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eBottom(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eRight(*(bc), &(*barrelRecHits), &(*topology)));
    bc_rook[bc_n] = *(max_element(rook_vect.begin(), rook_vect.end()));

    std::vector<float> vCov = EcalClusterTools::covariances( *(bc), &(*barrelRecHits), &(*topology), geometry);
    bc_see[bc_n] = vCov[0];
    bc_sep[bc_n] = vCov[1];
    bc_spp[bc_n] = vCov[2];


    // the indices in the following seem to be in the order: eta_min, eta_max, phi_min, phi_max
    // 
    // (where all numbers are integers and relative to )
    // 
    // see also http://cmslxr.fnal.gov/lxr/source/RecoEcal/EgammaCoreTools/src/EcalClusterTools.cc
    //
    //                                                                                                          eta        phi
    // 
    // 5 crystals in phi at deta = 0 (first line), -1 (second line), +1 (third line)
    bc_s1x5_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),     0,  0,    -2,  2);
    bc_s1x5_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -1, -1,    -2,  2);
    bc_s1x5_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),     1,  1,    -2,  2);
                                     
    // 3 crystals in phi at deta = 0 (first line), -1 (second line), +1 (third line)
    bc_s1x3_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),     0,  0,    -1,  1);
    bc_s1x3_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -1, -1,    -1,  1);
    bc_s1x3_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),     1,  1,    -1,  1);
                                                                                    
    // 5 crystals in eta at dphi = 0 (first line), -1 (second line), +1 (third line)
    bc_s5x1_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -2,  2,     0,  0);
    bc_s5x1_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -2,  2,    -1, -1);
    bc_s5x1_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -2,  2,     1,  1);
                                                                                    
    // 3 crystals in eta at dphi = 0 (first line), -1 (second line), +1 (third line)                                     
    bc_s3x1_0[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -1,  1,     0,  0);
    bc_s3x1_1[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -1,  1,    -1, -1);
    bc_s3x1_2[bc_n]=EcalClusterTools::matrixEnergy(*(bc), &(*barrelRecHits), &(*topology), mypair.first(),    -1,  1,     1,  1);

    // localCovariances no longer takes geometry as input.  At this point still need to do "cvs co RecoEcal/EgammaClusterTools" from the cern cvs. REMOVE COMENT ONCE IT WORKS
    bc_sieie[bc_n] = sqrt(EcalClusterTools::localCovariances(*(bc), &(*barrelRecHits), &(*topology))[0]);// WAS USING endcapRecHits (CAN SEE FROM LINE BELOW)
    //bc_sieie[bc_n] = sqrt(EcalClusterTools::localCovariances(*(bc), &(*endcapRecHits), &(*topology), &geometry)[0]);
    //bc_2x5_max[bc_n] = EcalClusterTools::e2x5Max(*(bc), &(*endcapRecHits), &(*topology));
    //bc_5x1_sam[bc_n] = EcalClusterTools::e5x1(*(bc), &(*endcapRecHits), &(*topology));
    bc_type[bc_n] = 1;

    bc_n++;
    bc_hybrid_n++;
  }
}
