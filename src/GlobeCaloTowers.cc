#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCaloTowers.h"


GlobeCaloTowers::GlobeCaloTowers(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  calotowerColl =  iConfig.getParameter<edm::InputTag>("CaloTowerColl");
  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeCaloTowers::defineBranch(TTree* tree) {
    
  tree->Branch("ct_n", &ct_n, "ct_n/I");
  
  ct_p4 = new TClonesArray("TLorentzVector", MAX_CALOTOWERS);
  tree->Branch("ct_p4", "TClonesArray", &ct_p4, 32000, 0);

  tree->Branch("ct_emEnergy", &ct_emEnergy, "ct_emEnergy[ct_n]/F");
  tree->Branch("ct_hadEnergy", &ct_hadEnergy, "ct_hadEnergy[ct_n]/F");
  tree->Branch("ct_outerEnergy", &ct_outerEnergy, "ct_outerEnergy[ct_n]/F");
  tree->Branch("ct_emL1", &ct_emL1, "ct_emL1[ct_n]/I");
  tree->Branch("ct_hadL1", &ct_hadL1, "ct_hadL1[ct_n]/I");
  tree->Branch("ct_size", &ct_size, "ct_size[ct_n]/I");
}

bool GlobeCaloTowers::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // get collections
  edm::Handle<edm::SortedCollection<CaloTower> > caltowH;
  iEvent.getByLabel(calotowerColl, caltowH);

  if (debug_level > 9)
    std::cout << "GlobeCaloTowers: CaloTower collection size: "<< caltowH->size() << std::endl;

  ct_p4->Clear();

  edm::ESHandle<CaloGeometry> pG;
  const CaloGeometry* caloGeom;
  iSetup.get<CaloGeometryRecord>().get(pG);
  caloGeom = pG.product();
  
  ct_n = 0;
  
  for(unsigned int i=0; i<caltowH->size(); i++) {

    if (ct_n >= MAX_CALOTOWERS) {
      std::cout << "GlobeCaloTowers: WARNING TOO MANY CALOTOWERS: " << caltowH->size() << " (allowed " << MAX_CALOTOWERS << ")" << std::endl;
      break;
    }

    CaloTowerRef ct(caltowH, i);
    // apply the cuts
	 if(gCUT->cut(*ct))continue; 
    // passed cuts

    GlobalPoint posi = caloGeom->getPosition(ct->id());
    new ((*ct_p4)[ct_n]) TLorentzVector();
    ((TLorentzVector *)ct_p4->At(ct_n))->SetXYZT(posi.x(), posi.y(), posi.z(), ct->energy());

    ct_emEnergy[ct_n] = ct->emEnergy();
    ct_hadEnergy[ct_n] = ct->hadEnergy();
    ct_outerEnergy[ct_n] = ct->outerEnergy();
    ct_emL1[ct_n] = ct->emLvl1();
    ct_hadL1[ct_n] = ct->hadLv11();
    ct_size[ct_n] = ct->constituentsSize();
    ct_n++;
  }

  return true;
}
