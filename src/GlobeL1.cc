#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeL1.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"

GlobeL1::GlobeL1(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  l1EMIso =  iConfig.getParameter<edm::InputTag>("L1EMIso");
  l1EMNonIso =  iConfig.getParameter<edm::InputTag>("L1EMNonIso");
  l1CenJet = iConfig.getParameter<edm::InputTag>("L1CenJet");
  l1ForJet = iConfig.getParameter<edm::InputTag>("L1ForJet");
  l1TauJet = iConfig.getParameter<edm::InputTag>("L1TauJet");
  l1EtMiss = iConfig.getParameter<edm::InputTag>("L1EtMiss");
  l1Mu = iConfig.getParameter<edm::InputTag>("L1Mu"); 
  m_gtReadoutRecord = iConfig.getParameter<edm::InputTag>("L1GtReadoutRecordTag");
  m_gtObjectMapRecord = iConfig.getParameter<edm::InputTag>("L1GtObjectMapRecordTag");

  debug_level = iConfig.getParameter<int>("Debug_Level");


  l1bits_phy = new std::vector<int>;
  l1bits_tec = new std::vector<int>;
  l1_labels = new std::map<std::string, int>;
}

void GlobeL1::defineBranch(TTree* tree) {
  
  tree->Branch("l1emiso_n", &l1emiso_n,"l1emiso_n/I");
  tree->Branch("l1emiso_eta", &l1emiso_eta,"l1emiso_eta[l1emiso_n]/F");
  tree->Branch("l1emiso_et", &l1emiso_et,"l1emiso_et[l1emiso_n]/F");
  tree->Branch("l1emiso_phi", &l1emiso_phi,"l1emiso_phi[l1emiso_n]/F");
  //tree->Branch("l1emiso_rank", &l1emiso_rank,"l1emiso_rank[l1emiso_n]/I");

  tree->Branch("l1emnoniso_n", &l1emnoniso_n,"l1emnoniso_n/I");
  tree->Branch("l1emnoniso_et", &l1emnoniso_et,"l1emnoniso_et[l1emnoniso_n]/F");
  tree->Branch("l1emnoniso_eta", &l1emnoniso_eta,"l1emnoniso_eta[l1emnoniso_n]/F");
  tree->Branch("l1emnoniso_phi", &l1emnoniso_phi,"l1emnoniso_phi[l1emnoniso_n]/F");
  //tree->Branch("l1emnoniso_rank", &l1emnoniso_rank,"l1emnoniso_rank[l1emnoniso_n]/I");

  tree->Branch("l1cenjet_n", &l1cenjet_n,"l1cenjet_n/I");
  tree->Branch("l1cenjet_et", &l1cenjet_et,"l1cenjet_et[l1cenjet_n]/F");
  tree->Branch("l1cenjet_eta", &l1cenjet_eta,"l1cenjet_eta[l1cenjet_n]/F");
  tree->Branch("l1cenjet_phi", &l1cenjet_phi,"l1cenjet_phi[l1cenjet_n]/F");
  //tree->Branch("l1cenjet_rank", &l1cenjet_rank,"l1cenjet_rank[l1cenjet_n]/I");

  tree->Branch("l1forjet_n", &l1forjet_n,"l1forjet_n/I");
  tree->Branch("l1forjet_et", &l1forjet_et,"l1forjet_et[l1forjet_n]/F");
  tree->Branch("l1forjet_eta", &l1forjet_eta,"l1forjet_eta[l1forjet_n]/F");
  tree->Branch("l1forjet_phi", &l1forjet_phi,"l1forjet_phi[l1forjet_n]/F");
  //tree->Branch("l1forjet_rank", &l1forjet_rank,"l1forjet_rank[l1forjet_n]/I");

  tree->Branch("l1taujet_n", &l1taujet_n,"l1taujet_n/I");
  tree->Branch("l1taujet_et", &l1taujet_et,"l1taujet_et[l1taujet_n]/F");
  tree->Branch("l1taujet_eta", &l1taujet_eta,"l1taujet_eta[l1taujet_n]/F");
  tree->Branch("l1taujet_phi", &l1taujet_phi,"l1taujet_phi[l1taujet_n]/F");
  //tree->Branch("l1taujet_rank", &l1taujet_rank,"l1taujet_rank[l1taujet_n]/I");

  //tree->Branch("l1met_n", &l1met_n,"l1met_n/I");
  tree->Branch("l1met_et", &l1met_et,"l1met_et/F");
  tree->Branch("l1met_phi", &l1met_phi,"l1met_phi/F");
  //tree->Branch("l1met_rank", &l1met_rank,"l1met_rank/I");

  tree->Branch("l1mu_n", &l1mu_n,"l1mu_n/I"); //CHECK, why sometimes there are 5? is it normal?
  tree->Branch("l1mu_et", &l1mu_et,"l1mu_et[l1mu_n]/F");
  tree->Branch("l1mu_eta", &l1mu_eta,"l1mu_eta[l1mu_n]/F");
  tree->Branch("l1mu_phi", &l1mu_phi,"l1mu_phi[l1mu_n]/F");
  //tree->Branch("l1mu_rank", &l1mu_rank,"l1mu_rank[l1mu_n]/I");
  
  tree->Branch("l1_labels", "std::map<std::string, int>", &l1_labels);
  tree->Branch("l1bits_phy", "std::vector<int>", &l1bits_phy);
  tree->Branch("l1bits_tec", "std::vector<int>", &l1bits_tec);
}

bool GlobeL1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  l1bits_phy->clear();
  l1bits_tec->clear();
  
   // access the L1 decisions
  edm::Handle<L1GlobalTriggerReadoutRecord> h_gtReadoutRecord;
  iEvent.getByLabel(m_gtReadoutRecord, h_gtReadoutRecord);
  
  edm::Handle<L1GlobalTriggerObjectMapRecord> h_gtObjectMapRecord;
  //iEvent.getByLabel(m_gtObjectMapRecord, h_gtObjectMapRecord);

  edm::Handle<l1extra::L1EmParticleCollection> l1emiso;
  edm::Handle<l1extra::L1EmParticleCollection> l1emnoniso;
  edm::Handle<l1extra::L1JetParticleCollection> l1cenjet;
  edm::Handle<l1extra::L1JetParticleCollection> l1forjet;
  edm::Handle<l1extra::L1JetParticleCollection> l1taujet;
  edm::Handle<l1extra::L1EtMissParticleCollection> l1met;
  edm::Handle<l1extra::L1MuonParticleCollection> l1mu;
  
  iEvent.getByLabel(l1EMIso, l1emiso);
  iEvent.getByLabel(l1EMNonIso, l1emnoniso);
  iEvent.getByLabel(l1CenJet, l1cenjet);
  iEvent.getByLabel(l1ForJet, l1forjet);
  iEvent.getByLabel(l1TauJet, l1taujet);
  iEvent.getByLabel(l1EtMiss, l1met);
  iEvent.getByLabel(l1Mu, l1mu);

  
  const std::vector<bool> & physics = h_gtReadoutRecord->decisionWord(0);
  for (unsigned int i = 0; i < 128; ++i) {
    l1bits_phy->push_back((int)physics[i]);
  }

  const std::vector<bool> & technical = h_gtReadoutRecord->technicalTriggerWord(0);
  for (unsigned int i = 0; i < 64; ++i) {
    l1bits_tec->push_back((int)technical[i]);
  }
  

  l1emiso_n = 0, l1emnoniso_n = 0, l1cenjet_n = 0, l1forjet_n = 0,
    l1taujet_n = 0, l1mu_n = 0;
  
  if (debug_level > 9){
    std::cout << "GlobeL1: L1EMIso collection size: "<< l1emiso->size() << std::endl;
    std::cout << "GlobeL1: L1EMNonIso collection size: "<< l1emnoniso->size() << std::endl;
    std::cout << "GlobeL1: L1CenJet collection size: "<< l1cenjet->size() << std::endl;
    std::cout << "GlobeL1: L1ForJet collection size: "<< l1forjet->size() << std::endl;
    std::cout << "GlobeL1: L1TauJet collection size: "<< l1taujet->size() << std::endl;
    std::cout << "GlobeL1: L1Mu collection size: "<< l1mu->size() << std::endl;
  }
  
  for(unsigned int i=0; i<l1emiso->size(); i++) {

    if (l1emiso_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1EMISO: " << l1emiso->size() << " (allowed " << MAX_L1 << ")" << std::endl;
      break;
    }

    l1extra::L1EmParticleRef iso(l1emiso, i);
    l1emiso_eta[i] = iso->eta();
    l1emiso_phi[i] = iso->phi();
    l1emiso_et[i] = iso->et();
    //      l1emiso_rank[i] = iso->rank(); // 
      
    l1emiso_n++;
  }

  for(unsigned int i=0; i<l1emnoniso->size(); i++) {
      
    if (l1emnoniso_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1EMNONISO: " << l1emnoniso->size() << " (allowed " << MAX_L1 << ")" << std::endl;
      break;
    }

    l1extra::L1EmParticleRef noniso(l1emnoniso, i);
      
    l1emnoniso_eta[i] = noniso->eta();
    l1emnoniso_phi[i] = noniso->phi();
    l1emnoniso_et[i] = noniso->et();
    l1emnoniso_n++;
  }

  for(unsigned int i=0; i<l1cenjet->size(); i++) {
    
    if (l1cenjet_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1CENJET: " << l1cenjet->size() << " (allowed " << MAX_L1 << ")" << std::endl;
      break;
    }
    l1extra::L1JetParticleRef cenjet(l1cenjet, i);
    l1cenjet_eta[i] = cenjet->eta();
    l1cenjet_phi[i] = cenjet->phi();
    l1cenjet_et[i] = cenjet->et();
    l1cenjet_n++;
  }
 
  for(unsigned int i=0; i<l1forjet->size(); i++) {
    
    if (l1forjet_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1FORJET: " << l1forjet->size() << " (allowed " << MAX_L1<< ")" << std::endl;
      break;
    }
    l1extra::L1JetParticleRef forjet(l1forjet, i);
    l1forjet_eta[i] = forjet->eta();
    l1forjet_phi[i] = forjet->phi();
    l1forjet_et[i] = forjet->et();
    l1forjet_n++;
  }

  for(unsigned int i=0; i<l1taujet->size(); i++) {
    
    if (l1taujet_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1TAUJET: " << l1taujet->size() << " (allowed " << MAX_L1 << ")" << std::endl;
      break;
    }

    l1extra::L1JetParticleRef taujet(l1taujet, i);
    l1taujet_eta[i] = taujet->eta();
    l1taujet_phi[i] = taujet->phi();
    l1taujet_et[i] = taujet->et();
    l1taujet_n++;
  }

  for(unsigned int i=0; i<l1mu->size(); i++) {
    
    if (l1mu_n >= MAX_L1) {
      std::cout << "GlobeL1: WARNING TOO MANY L1MU: " << l1mu->size() << " (allowed " << MAX_L1 << ")" << std::endl;
      break;
    }

    l1extra::L1MuonParticleRef mu(l1mu, i);
    l1mu_eta[i] = mu->eta();
    l1mu_phi[i] = mu->phi();
    l1mu_et[i] = mu->et();
    l1mu_n++;
  }

  if (l1met->size() > 0) {
    l1met_et = l1met->begin()->et();
    l1met_phi = l1met->begin()->phi();
  } else {
    l1met_et = 0;
    l1met_phi = 0;
  }
 
  return true;
}
