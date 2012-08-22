#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalHits.h"

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
#include <numeric>

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

  CrackCorrFunc = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
  LocalCorrFunc = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);
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
  tree->Branch("sc_sphi", &sc_sphi, "sc_sphi[sc_n]/F");
  tree->Branch("sc_seta", &sc_seta, "sc_seta[sc_n]/F");
  tree->Branch("sc_brem", &sc_brem, "sc_brem[sc_n]/F");
  tree->Branch("sc_r9", &sc_r9, "sc_r9[sc_n]/F");
  tree->Branch("sc_nbc", &sc_nbc, "sc_nbc[sc_n]/I");
  tree->Branch("sc_bcseedind", &sc_bcseedind, "sc_bcseedind[sc_n]/I");
  sprintf (a2, "sc_bcind[sc_n][%d]/I", MAX_SUPERCLUSTER_BASICCLUSTERS);
  tree->Branch("sc_bcind", &sc_bcind, a2);
  sprintf (a2, "sc_bccrackcorr[sc_n][%d]/F", MAX_SUPERCLUSTER_BASICCLUSTERS);
  tree->Branch("sc_bccrackcorr",&sc_bccrackcorr, a2);
  sprintf (a2, "sc_bclocalcorr[sc_n][%d]/F", MAX_SUPERCLUSTER_BASICCLUSTERS);
  tree->Branch("sc_bclocalcorr",&sc_bclocalcorr, a2);

  //basic clusters
  tree->Branch("bc_n", &bc_n, "bc_n/I");
  tree->Branch("bc_islbar_n", &bc_islbar_n, "bc_islbar_n/I");
  tree->Branch("bc_islend_n", &bc_islend_n, "bc_islend_n/I");
  tree->Branch("bc_hybrid_n", &bc_hybrid_n, "bc_hybrid_n/I");
  tree->Branch("bc_p4", "TClonesArray", &bc_p4, 32000, 0);
  tree->Branch("bc_xyz", "TClonesArray", &bc_xyz, 32000, 0);
  tree->Branch("bc_nhits", &bc_nhits,"bc_nhits[bc_n]/I");
  tree->Branch("bc_s1", &bc_s1, "bc_s1[bc_n]/F");
  tree->Branch("bc_chx", &bc_chx, "bc_chx[bc_n]/F");
  tree->Branch("bc_s4", &bc_s4, "bc_s4[bc_n]/F");
  tree->Branch("bc_s9", &bc_s9, "bc_s9[bc_n]/F");
  tree->Branch("bc_s25", &bc_s25, "bc_s25[bc_n]/F");
  tree->Branch("bc_sipip", &bc_sipip, "bc_sipip[bc_n]/F");
  tree->Branch("bc_sieie", &bc_sieie, "bc_sieie[bc_n]/F");
  tree->Branch("bc_sieip", &bc_sieip, "bc_sieip[bc_n]/F");
  tree->Branch("bc_type", &bc_type, "bc_type[bc_n]/I");//type 1 = hybrid, 2 = island endcap, 3 = island barrel.
  //tree->Branch("bc_seed", &bc_seed, "bc_seed[bc_n]/I");
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

  //const std::string instan="hybridBarrelBasicClusters";
  iEvent.getByLabel(barrelHybridClusterColl, hybridClustersBarrelH);
  if (barrelBasicClusterColl.encode() != "") 
    iEvent.getByLabel(barrelBasicClusterColl, basicClustersBarrelH);

  iEvent.getByLabel(endcapBasicClusterColl, basicClustersEndcapH);

  CrackCorrFunc->init(iSetup);
  LocalCorrFunc->init(iSetup);

  if (debug_level > 9) {
    std::cout << "GlobeEcalClusters: superClustersEndcapH->size() "<< superClustersEndcapH->size() << std::endl;
    std::cout << "GlobeEcalClusters: hybridClustersBarrelH->size() "<< hybridClustersBarrelH->size() << std::endl;
    std::cout << "GlobeEcalClusters: basicClustersEndcapH->size() "<< basicClustersEndcapH->size() << std::endl;
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
    
  return true;
}

//----------------------------------------------------------------------

void GlobeEcalClusters::e2xNOe5xN(float& e2xN, float& e5xN, const reco::SuperCluster* superCluster, const EcalRecHitCollection* hits,
                                  const  CaloSubdetectorTopology* topology) {
  
  std::vector<std::pair<DetId, float> > crystals = superCluster->hitsAndFractions();

  // look for the max energy crystal
  double eMax=0;
  DetId eMaxId(0);
  std::vector<std::pair<DetId, float> >::iterator it;
  EcalRecHitCollection::const_iterator itt;
  
  EcalRecHit testEcalRecHit;
  
  for(it = crystals.begin(); it != crystals.end(); it++) {
    
    itt = hits->find((*it).first);
    if (((*it).first != DetId(0)) && (itt != hits->end())) {
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
GlobeEcalClusters::analyzeBarrelSuperClusters() {

  for(unsigned int i=0; i<superClustersHybridH->size(); i++) {
    
    if(sc_n >= MAX_SUPERCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many superclusters. " << MAX_SUPERCLUSTERS << " allowed.  Event has " << superClustersHybridH->size() +  superClustersEndcapH->size()<< std::endl;
      break;
    }
    
    reco::SuperClusterRef sc(superClustersHybridH, i);

    sc_bcseedind[sc_n] = -1;
    // get index to seed basic cluster
    for(unsigned int j=0; j<hybridClustersBarrelH->size(); ++j) {
      reco::BasicClusterRef basic(hybridClustersBarrelH, j);
      if (&(*sc->seed()) == &(*basic)) {
        //sc_bcseedind[sc_n] = j + basicClustersEndcapH->size();
	sc_bcseedind[sc_n] = j + bc_islend_n;
        break;
      }
    }
    
    if (sc_bcseedind[sc_n] == -1)
      continue;

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
    sc_sieie[sc_n] = sqrt(EcalClusterTools::scLocalCovariances(*(sc), &(*barrelRecHits), &(*topology))[0]);

    // sigma_phi and sigma_eta as defined in reco::SuperCluster
    sc_sphi[sc_n] = sc->phiWidth();
    sc_seta[sc_n] = sc->etaWidth();
    if (sc_seta[sc_n]>0) sc_brem[sc_n] = sc_sphi[sc_n] / sc_seta[sc_n];
    else sc_brem[sc_n]=-1; // something is not ok with sigma_eta, in this case

    // SC r9
    if (sc->rawEnergy()>0) sc_r9[sc_n] = EcalClusterTools::e3x3(  *(sc->seed()), &(*barrelRecHits), &(*topology)) / sc->rawEnergy();
    else sc_r9[sc_n]=-1;

    //SEED BC 
    if (debug_level > 10)
      std::cout << sc->clustersSize() << std::endl;
    sc_nbc[sc_n]= sc->clustersSize();
    
    // get indices to basic clusters
    if (sc->clustersSize() > 0) { 
      int limit = 0;
      
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
        
        if (limit >= MAX_SUPERCLUSTER_BASICCLUSTERS) {
          std::cout << "GlobeEcalCluster: WARNING too many basiclusters. in basicClustersBarrelH (" << 
            MAX_SUPERCLUSTER_BASICCLUSTERS << " allowed). Event has " << sc->clustersSize() << std::endl;
          break;
        } 
        
        for(unsigned int j=0; j<hybridClustersBarrelH->size(); ++j) {
          reco::BasicClusterRef basic(hybridClustersBarrelH, j);
          if (&(**itClus) == &(*basic)) {
            //sc_bcind[sc_n][limit] = j + basicClustersEndcapH->size();
	    sc_bcind[sc_n][limit] = j + bc_islend_n;
            break;
          }
        }
	
	const reco::CaloClusterPtr cc = *itClus;
	sc_bccrackcorr[sc_n][limit] = CrackCorrFunc->getValue(*cc);
	sc_bclocalcorr[sc_n][limit] = LocalCorrFunc->getValue(*cc);

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
 
    // sigma_phi and sigma_eta as defined in reco::SuperCluster
    sc_sphi[sc_n] = sc->phiWidth();
    sc_seta[sc_n] = sc->etaWidth();
    if (sc_seta[sc_n]>0) sc_brem[sc_n] = sc_sphi[sc_n] / sc_seta[sc_n];
    else sc_brem[sc_n]=-1; // something is not ok with sigma_eta, in this case

    // SC r9
    if (sc->rawEnergy()>0) sc_r9[sc_n] = EcalClusterTools::e3x3(  *(sc->seed()), &(*endcapRecHits), &(*topology)) / sc->rawEnergy();
    else sc_r9[sc_n]=-1;

    // get index to seed basic cluster
    for(unsigned int j=0; j<basicClustersEndcapH->size(); ++j) {
      reco::BasicClusterRef basic(basicClustersEndcapH, j);

       if(gCUT->cut(*basic))
	 continue;
      if (&(*sc->seed()) == &(*basic)) {
        sc_bcseedind[sc_n] = j;
        break;
      }
    }
    
    // get indices to basic clusters
    if (sc->clustersSize() > 0) { 
      int limit = 0;
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
        if (limit >= MAX_SUPERCLUSTER_BASICCLUSTERS) {
          std::cout << "GlobeEcalCluster: WARNING too many basiclusters. in basicClustersEndcapH (" << MAX_SUPERCLUSTER_BASICCLUSTERS << " allowed). Event has " << sc->clustersSize() << std::endl;
          break;
        } 
        for(unsigned int j=0; j<basicClustersEndcapH->size(); ++j) {
          reco::BasicClusterRef basic(basicClustersEndcapH, j);
	  if(gCUT->cut(*basic))
	    continue;
	  if (&(**itClus) == &(*basic)) {
            sc_bcind[sc_n][limit] = j;
            break;
          }
        }

	const reco::CaloClusterPtr cc = *itClus;
	sc_bccrackcorr[sc_n][limit] = CrackCorrFunc->getValue(*cc);
	sc_bclocalcorr[sc_n][limit] = LocalCorrFunc->getValue(*cc);

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
GlobeEcalClusters::analyzeEndcapBasicClusters() {

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

    new ((*bc_xyz)[bc_n]) TVector3();
    ((TVector3 *)bc_xyz->At(bc_n))->SetXYZ(bc->position().x(), bc->position().y(), bc->position().z());

    std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
    bc_nhits[bc_n] = hits.size(); 

    bc_seed[bc_n] = -1;
    for(EERecHitCollection::const_iterator it=endcapRecHits->begin(); it!=endcapRecHits->end(); it++) {
      if (bc->seed().rawId() == it->detid().rawId())
        bc_seed[bc_n] = (it - endcapRecHits->begin()); 
    }

    bc_s1[bc_n] = EcalClusterTools::eMax(*(bc), &(*endcapRecHits)); 
    bc_s4[bc_n] = EcalClusterTools::e2x2(*(bc), &(*endcapRecHits), &(*topology)); 
    bc_s9[bc_n] = EcalClusterTools::e3x3(*(bc), &(*endcapRecHits), &(*topology)); 
    bc_s25[bc_n] = EcalClusterTools::e5x5(*(bc), &(*endcapRecHits), &(*topology)); 
    
    std::vector<float> rook_vect;
    rook_vect.push_back(EcalClusterTools::eLeft(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eTop(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eBottom(*(bc), &(*endcapRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eRight(*(bc), &(*endcapRecHits), &(*topology)));
    bc_chx[bc_n] = std::accumulate(rook_vect.begin(), rook_vect.end(), 0.);

    std::vector<float> vCov = EcalClusterTools::localCovariances( *(bc), &(*endcapRecHits), &(*topology));
    bc_sieie[bc_n] = sqrt(vCov[0]);
    bc_sieip[bc_n] = vCov[1];
    bc_sipip[bc_n] = sqrt(vCov[2]);

    bc_type[bc_n] = 2;
    bc_n++;
    bc_islend_n++;
  }
}

//----------------------------------------------------------------------
 void GlobeEcalClusters::analyzeBarrelHybridClusters() { 

   if (barrelBasicClusterColl.encode() != "") {

     // BASIC CLUSTERS BARREL //CHECK: missing type of clusters //MARCO
     for(unsigned int i=0; i< basicClustersBarrelH->size(); i++) {
       
       if(bc_n >= MAX_BASICCLUSTERS) {
         std::cout << "GlobeEcalCluster: WARNING too many basicclusters. " << MAX_BASICCLUSTERS << " allowed.  Event has " << basicClustersEndcapH->size() + hybridClustersBarrelH->size() + basicClustersBarrelH->size()<< std::endl;
         break;
       }
       
       reco::BasicClusterRef bc(basicClustersBarrelH, i);
       
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
       
       new ((*bc_xyz)[bc_n]) TVector3();
       ((TVector3 *)bc_xyz->At(bc_n))->SetXYZ(bc->position().x(), bc->position().y(), bc->position().z());
       
       std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
       bc_nhits[bc_n] = hits.size(); 
       
       bc_seed[bc_n] = -1;
       for(EBRecHitCollection::const_iterator it=barrelRecHits->begin(); it!=barrelRecHits->end(); it++) {
         if (bc->seed().rawId() == it->detid().rawId())
           bc_seed[bc_n] = (it - barrelRecHits->begin()); 
       }

       bc_s1[bc_n] = EcalClusterTools::eMax(*(bc), &(*barrelRecHits)); 
       bc_s4[bc_n] = EcalClusterTools::e2x2(*(bc), &(*barrelRecHits), &(*topology)); 
       bc_s9[bc_n] = EcalClusterTools::e3x3(*(bc), &(*barrelRecHits), &(*topology)); 
       bc_s25[bc_n] = EcalClusterTools::e5x5(*(bc), &(*barrelRecHits), &(*topology)); 
       
       std::vector<float> rook_vect;
       rook_vect.push_back(EcalClusterTools::eLeft(*(bc), &(*barrelRecHits), &(*topology)));
       rook_vect.push_back(EcalClusterTools::eTop(*(bc), &(*barrelRecHits), &(*topology)));
       rook_vect.push_back(EcalClusterTools::eBottom(*(bc), &(*barrelRecHits), &(*topology)));
       rook_vect.push_back(EcalClusterTools::eRight(*(bc), &(*barrelRecHits), &(*topology)));
       //bc_rook[bc_n] = *(max_element(rook_vect.begin(), rook_vect.end()));
       bc_chx[bc_n] = std::accumulate(rook_vect.begin(), rook_vect.end(), 0.);
       
       std::vector<float> vCov = EcalClusterTools::localCovariances( *(bc), &(*barrelRecHits), &(*topology));
       
       bc_sieie[bc_n] = sqrt(vCov[0]);
       bc_sieip[bc_n] = vCov[1];
       bc_sipip[bc_n] = sqrt(vCov[2]);

       bc_type[bc_n] = 3;
       bc_n++;
       bc_islbar_n++;
     }
   }
 
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
       
    new ((*bc_xyz)[bc_n]) TVector3();
    ((TVector3 *)bc_xyz->At(bc_n))->SetXYZ(bc->position().x(), bc->position().y(), bc->position().z());

    std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
    bc_nhits[bc_n] = hits.size(); // CHECK no direct method in the dataFormat only getRecHitsByDetId

    bc_seed[bc_n] = -1;
    for(EBRecHitCollection::const_iterator it=barrelRecHits->begin(); it!=barrelRecHits->end(); it++) {
      if (bc->seed().rawId() == it->detid().rawId())
        bc_seed[bc_n] = (it - barrelRecHits->begin()); 
    }
    
    //std::pair<DetId, float> mypair=EcalClusterTools::getMaximum( *(bc), &(*barrelRecHits)); 
    bc_s1[bc_n] = EcalClusterTools::eMax(  *(bc), &(*barrelRecHits));
    bc_s4[bc_n] = EcalClusterTools::e2x2(  *(bc), &(*barrelRecHits), &(*topology));
    bc_s9[bc_n] = EcalClusterTools::e3x3(  *(bc), &(*barrelRecHits), &(*topology));
    bc_s25[bc_n] = EcalClusterTools::e5x5(  *(bc), &(*barrelRecHits), &(*topology));

    std::vector<float> rook_vect;
    rook_vect.push_back(EcalClusterTools::eLeft(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eTop(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eBottom(*(bc), &(*barrelRecHits), &(*topology)));
    rook_vect.push_back(EcalClusterTools::eRight(*(bc), &(*barrelRecHits), &(*topology)));
    //bc_rook[bc_n] = *(max_element(rook_vect.begin(), rook_vect.end()));
    bc_chx[bc_n] = std::accumulate(rook_vect.begin(), rook_vect.end(), 0.);

    std::vector<float> vCov = EcalClusterTools::localCovariances( *(bc), &(*barrelRecHits), &(*topology));
    bc_sieie[bc_n] = sqrt(vCov[0]);
    bc_sieip[bc_n] = vCov[1];
    bc_sipip[bc_n] = sqrt(vCov[2]);

    //bc_sieie[bc_n] = sqrt(EcalClusterTools::localCovariances(*(bc), &(*endcapRecHits), &(*topology), &geometry)[0]);
    //bc_2x5_max[bc_n] = EcalClusterTools::e2x5Max(*(bc), &(*endcapRecHits), &(*topology));
    //bc_5x1_sam[bc_n] = EcalClusterTools::e5x1(*(bc), &(*endcapRecHits), &(*topology));
    
    bc_type[bc_n] = 1;
    bc_n++;
    bc_hybrid_n++;
  }
}

std::vector<float> GlobeEcalClusters::getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row) {
  std::vector<float> esHits;

  const GlobalPoint point(X,Y,Z);

  const CaloSubdetectorGeometry *geometry_p ;
  geometry_p = geometry.getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  std::map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;
    
  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
    
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

std::vector<float> GlobeEcalClusters::getESShape(std::vector<float> ESHits0)
{
  std::vector<float> esShape;
    
  const int nBIN = 21;
  float esRH_F[nBIN];
  float esRH_R[nBIN];
  for (int idx=0; idx<nBIN; idx++) {
    esRH_F[idx] = 0.;
    esRH_R[idx] = 0.;
  }

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      esRH_F[(nBIN-1)/2] = ESHits0[ibin];
      esRH_R[(nBIN-1)/2] = ESHits0[ibin+31];
    } else {
      esRH_F[(nBIN-1)/2+ibin] = ESHits0[ibin];
      esRH_F[(nBIN-1)/2-ibin] = ESHits0[ibin+15];
      esRH_R[(nBIN-1)/2+ibin] = ESHits0[ibin+31];
      esRH_R[(nBIN-1)/2-ibin] = ESHits0[ibin+31+15];
    }
  } 

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=0; id_X<21; id_X++) {
    totalEnergyXX  += esRH_F[id_X];
    EffStatsXX     += esRH_F[id_X]*(id_X-10)*(id_X-10);
    totalEnergyYY  += esRH_R[id_X];
    EffStatsYY     += esRH_R[id_X]*(id_X-10)*(id_X-10);
  }
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;

  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);
    
  return esShape;
}
