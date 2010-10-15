#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHLT.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TLorentzVector.h"

#include <bitset>

void set(unsigned int& x, int bit) {
  x |= (1 << bit);
}

void reset(unsigned int& x, int bit) {
  x &= ~ (1 << bit);
}

bool check(unsigned int x, int bit) {
  return (x && (1 << bit));
}

GlobeHLT::GlobeHLT(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  edm::ParameterSet psetHLT = iConfig.getParameter<edm::ParameterSet>("HLTParameters");
  fullHLT   = psetHLT.getParameter<bool>("FullHLT");
  inputTag_ = psetHLT.getParameter<edm::InputTag>("TriggerResultsTag");
  secondaryTriggerON  = psetHLT.getParameter<bool>("useSecondaryTrigger");
  
  // fullHLT
  if(fullHLT){
    theElHLTLabels = iConfig.getParameter<std::vector<edm::InputTag> >("ElectronHLTLabels");
    thePhHLTLabels = iConfig.getParameter<std::vector<edm::InputTag> >("PhotonHLTLabels");
    theMuHLTLabels = iConfig.getParameter<std::vector<edm::InputTag> >("MuonHLTLabels");
  } else {
    // new
    hlt1Tag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");
    hlt_path_names_HLT1_1 = new std::vector<std::string>; hlt_path_names_HLT1_1->clear();
    hlt_path_names_HLT1_2 = new std::vector<std::string>; hlt_path_names_HLT1_2->clear();
    hlt_path_names_HLT1_3 = new std::vector<std::string>; hlt_path_names_HLT1_3->clear();
    hlt_path_names_HLT1_4 = new std::vector<std::string>; hlt_path_names_HLT1_4->clear();
    if(secondaryTriggerON){
      hlt2Tag_  = psetHLT.getParameter<edm::InputTag>("SecondaryTriggerResultsTag");
      hlt_path_names_HLT2_1 = new std::vector<std::string>; hlt_path_names_HLT2_1->clear();
      hlt_path_names_HLT2_2 = new std::vector<std::string>; hlt_path_names_HLT2_2->clear();
      hlt_path_names_HLT2_3 = new std::vector<std::string>; hlt_path_names_HLT2_3->clear();
      hlt_path_names_HLT2_4 = new std::vector<std::string>; hlt_path_names_HLT2_4->clear();
    }
    theHLTLabels.clear();
    hlt_label_names_1 = new std::vector<std::string>; hlt_label_names_1->clear();
    hlt_label_names_2 = new std::vector<std::string>; hlt_label_names_2->clear();
    hlt_label_names_3 = new std::vector<std::string>; hlt_label_names_3->clear();
    hlt_label_names_4 = new std::vector<std::string>; hlt_label_names_4->clear();
  }
  
  electronColl   = iConfig.getParameter<edm::InputTag>("ElectronColl_std");
  photonCollStd  = iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  muonColl       = iConfig.getParameter<edm::InputTag>("MuonColl");
  jetColl        = iConfig.getParameter<edm::InputTag>("JetColl_it5");
  debug_level    = iConfig.getParameter<int>("Debug_Level");
  
  gCUT           = new GlobeCuts(iConfig);
}

void GlobeHLT::defineBranch(TTree* tree) {
  
  hlt_p4  = new TClonesArray("TLorentzVector",   MAX_HLT);
  
  if (fullHLT) {
    hlt_el_p4 =   new TClonesArray("TLorentzVector", MAX_HLT);
    hlt_ph_p4 =   new TClonesArray("TLorentzVector", MAX_HLT);
    hlt_mu_p4 =   new TClonesArray("TLorentzVector", MAX_HLT);
    hlt_el_tkp4 = new TClonesArray("TLorentzVector", MAX_HLT);
    tree->Branch("hlt_el_tkp4", "TClonesArray", &hlt_el_tkp4, 32000, 0);
    tree->Branch("hlt_el_n", &hlt_el_n, "hlt_el_n/I");
    tree->Branch("hlt_el_p4", "TClonesArray", &hlt_el_p4, 32000, 0);
    tree->Branch("hlt_el_offlineind", &hlt_el_offlineind, "hlt_el_offlineind[hlt_el_n]/I");
    tree->Branch("hlt_el_candpath", &hlt_el_candpath, "hlt_el_candpath[hlt_el_n]/i");
    
    tree->Branch("hlt_mu_n", &hlt_mu_n, "hlt_mu_n/I");
    tree->Branch("hlt_mu_p4", "TClonesArray", &hlt_mu_p4, 32000, 0);
    tree->Branch("hlt_mu_offlineind", &hlt_mu_offlineind, "hlt_mu_offlineind[hlt_mu_n]/I");
    tree->Branch("hlt_mu_candpath", &hlt_mu_candpath, "hlt_mu_candpath[hlt_mu_n]/i");
    
    tree->Branch("hlt_ph_n", &hlt_ph_n, "hlt_ph_n/I");
    tree->Branch("hlt_ph_p4", "TClonesArray", &hlt_ph_p4, 32000, 0);
    tree->Branch("hlt_ph_offlineind", &hlt_ph_offlineind, "hlt_ph_offlineind[hlt_ph_n]/I");
    tree->Branch("hlt_ph_candpath", &hlt_ph_candpath, "hlt_ph_candpath[hlt_ph_n]/i");
  } else {
    // Event Trigger
    tree->Branch("hlt1_bit_1", &hlt1_bit_1, "hlt1_bit_1/i");
    tree->Branch("hlt1_bit_2", &hlt1_bit_2, "hlt1_bit_2/i");
    tree->Branch("hlt1_bit_3", &hlt1_bit_3, "hlt1_bit_3/i");
    tree->Branch("hlt1_bit_4", &hlt1_bit_4, "hlt1_bit_4/i");
    tree->Branch("hlt_path_names_HLT1_1", "std::vector<std::string>", &hlt_path_names_HLT1_1);
    tree->Branch("hlt_path_names_HLT1_2", "std::vector<std::string>", &hlt_path_names_HLT1_2);
    tree->Branch("hlt_path_names_HLT1_3", "std::vector<std::string>", &hlt_path_names_HLT1_3);
    tree->Branch("hlt_path_names_HLT1_4", "std::vector<std::string>", &hlt_path_names_HLT1_4);
    //
    if(secondaryTriggerON){
      tree->Branch("hlt2_bit_1", &hlt2_bit_1, "hlt2_bit_1/i");
      tree->Branch("hlt2_bit_2", &hlt2_bit_2, "hlt2_bit_2/i");
      tree->Branch("hlt2_bit_3", &hlt2_bit_3, "hlt2_bit_3/i");
      tree->Branch("hlt2_bit_4", &hlt2_bit_4, "hlt2_bit_4/i");
      tree->Branch("hlt_path_names_HLT2_1", "std::vector<std::string>", &hlt_path_names_HLT2_1);
      tree->Branch("hlt_path_names_HLT2_2", "std::vector<std::string>", &hlt_path_names_HLT2_2);
      tree->Branch("hlt_path_names_HLT2_3", "std::vector<std::string>", &hlt_path_names_HLT2_3);
      tree->Branch("hlt_path_names_HLT2_4", "std::vector<std::string>", &hlt_path_names_HLT2_4);
    }
    // Trigger Candidates
    tree->Branch("hlt_n", &hlt_n, "hlt_n/I");
    tree->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0);
    tree->Branch("hlt_candpath_1"    , &hlt_candpath_1    , "hlt_candpath_1[hlt_n]/i"    );
    tree->Branch("hlt_candpath_2"    , &hlt_candpath_2    , "hlt_candpath_2[hlt_n]/i"    );
    tree->Branch("hlt_candpath_3"    , &hlt_candpath_3    , "hlt_candpath_3[hlt_n]/i"    );
    tree->Branch("hlt_candpath_4"    , &hlt_candpath_4    , "hlt_candpath_4[hlt_n]/i"    );
    tree->Branch("hlt_id"            , &hlt_id            , "hlt_id[hlt_n]/I"            );
    tree->Branch("hlt_el_offlineind" , &hlt_el_offlineind , "hlt_el_offlineind[hlt_n]/I" );
    tree->Branch("hlt_ph_offlineind" , &hlt_ph_offlineind , "hlt_ph_offlineind[hlt_n]/I" );
    tree->Branch("hlt_mu_offlineind" , &hlt_mu_offlineind , "hlt_mu_offlineind[hlt_n]/I" );
    tree->Branch("hlt_jet_offlineind", &hlt_jet_offlineind, "hlt_jet_offlineind[hlt_n]/I");
    //
    tree->Branch("hlt_label_names_1", "std::vector<std::string>", &hlt_label_names_1);
    tree->Branch("hlt_label_names_2", "std::vector<std::string>", &hlt_label_names_2);
    tree->Branch("hlt_label_names_3", "std::vector<std::string>", &hlt_label_names_3);
    tree->Branch("hlt_label_names_4", "std::vector<std::string>", &hlt_label_names_4);
  }
  
}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(electronColl, elH);
  
  edm::Handle<reco::MuonCollection> muH;
  iEvent.getByLabel(muonColl, muH);
  
  edm::Handle<reco::PhotonCollection> phH;
  iEvent.getByLabel(photonCollStd, phH);
  
  edm::Handle<reco::CaloJetCollection> jetH;
  iEvent.getByLabel(jetColl, jetH);  
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  edm::Handle<trigger::TriggerEventWithRefs> triggerObjWithRef;
  
  if (fullHLT)
    iEvent.getByLabel(inputTag_, triggerObjWithRef);
  else    
    iEvent.getByLabel(inputTag_, triggerObj);
  
  
  if (fullHLT) {
    hlt_el_p4->Clear();
    hlt_el_tkp4->Clear();
    hlt_el_n = 0;
    hlt_ph_p4->Clear();
    hlt_ph_n = 0;
    hlt_mu_p4->Clear();
    hlt_mu_n = 0;
    for(unsigned int n=0; n<theElHLTLabels.size(); n++) {
      std::vector<edm::Ref<reco::ElectronCollection> > recoecalcands;
      
      trigger::size_type elIndex = triggerObjWithRef->filterIndex(theElHLTLabels[n]);
      
      if (!(elIndex >= triggerObjWithRef->size())) {
	
	triggerObjWithRef->getObjects(elIndex, trigger::TriggerElectron, recoecalcands);
	
	for (unsigned int i=0; i<recoecalcands.size(); i++) {
	  
	  if (hlt_el_n >= MAX_HLT) {
	    std::cout << "GlobeHLT: WARNING TOO MANY ELECTRON HLT CANDIDATE: (allowed " << MAX_HLT << ")" << std::endl;
	    break;
	  }
	  
	  int isStored = -1;
	  TLorentzVector lv(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->energy());
	  for(int j=0; j<hlt_el_n; j++) {
	    TLorentzVector* lv_temp =  (TLorentzVector*)hlt_el_p4->At(j);
	    if (lv == *lv_temp) {
	      isStored = j;
	      break;
	    }
	  }
	  
	  if (isStored == -1) {
	    new ((*hlt_el_p4)[hlt_el_n]) TLorentzVector();
	    ((TLorentzVector *)hlt_el_p4->At(hlt_el_n))->SetXYZT(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->energy());
	    
	    new ((*hlt_el_tkp4)[hlt_el_n]) TLorentzVector();
	    ((TLorentzVector *)hlt_el_tkp4->At(hlt_el_n))->SetXYZT(recoecalcands[i]->track()->px(), recoecalcands[i]->track()->py(), recoecalcands[i]->track()->pz(), recoecalcands[i]->track()->p());
	    
	    float pt = recoecalcands[i]->pt();
	    float eta = recoecalcands[i]->eta();
	    float phi = recoecalcands[i]->phi();
	    TVector3 v1;
	    v1.SetPtEtaPhi(pt, eta, phi);
	    
	    float dRmin = 1.0;
	    int index = -1;
	    for(unsigned int z=0; z<elH->size(); z++) {
	      reco::GsfElectronRef e(elH, z);
	      // apply the cuts
	      if(gCUT->cut(*e))
		continue;
	      // passed cuts
	      
	      TVector3 v2;
	      v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	      float dR = v1.DeltaR(v2);
	      if (dR < dRmin) {
		dRmin = dR;
		index = z;
	      }
	    }
	    
	    hlt_el_offlineind[hlt_el_n] = index;
	    hlt_el_candpath[hlt_el_n] = 0;
	    set(hlt_el_candpath[hlt_el_n], n);
	    hlt_el_n++;
	  } else {
	    set(hlt_el_candpath[isStored], n);
	  }
	}
      }
    }
    
    // PHOTON TRIGGER
    for(unsigned int n=0; n<thePhHLTLabels.size(); n++) {
      std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > recoecalcands;
      trigger::size_type elIndex = triggerObjWithRef->filterIndex(thePhHLTLabels[n]);
      
      if (!(elIndex >= triggerObjWithRef->size())) {
	triggerObjWithRef->getObjects(elIndex, trigger::TriggerPhoton, recoecalcands);
	
	for (unsigned int i=0; i<recoecalcands.size(); i++) {
	  
	  if (hlt_ph_n >= MAX_HLT) {
	    std::cout << "GlobeHLT: WARNING TOO MANY PHOTON HLT CANDIDATE: (allowed " << MAX_HLT << ")" << std::endl;
	    break;
	  }
	  
	  int isStored = -1;
	  TLorentzVector lv(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->energy());
	  for(int j=0; j<hlt_ph_n; j++) {
	    TLorentzVector* lv_temp =  (TLorentzVector*)hlt_ph_p4->At(j);
	    if (lv == *lv_temp) {
	      isStored = j;
	      break;
	    }
	  }
	  
	  if (isStored == -1) {
	    new ((*hlt_ph_p4)[hlt_ph_n]) TLorentzVector();
	    ((TLorentzVector *)hlt_ph_p4->At(hlt_ph_n))->SetXYZT(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->energy());
	    
	    float pt = recoecalcands[i]->pt();
	    float eta = recoecalcands[i]->eta();
	    float phi = recoecalcands[i]->phi();
	    TVector3 v1;
	    v1.SetPtEtaPhi(pt, eta, phi);
	    
	    float dRmin = 1.0;
	    int index = -1;
	    for(unsigned int z=0; z<phH->size(); z++) {
	      reco::PhotonRef e(phH, z);
	      // apply the cuts
	      if(gCUT->cut(*e))
		continue;
	      // passed cuts
	      
	      TVector3 v2;
	      v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	      float dR = v1.DeltaR(v2);
	      if (dR < dRmin) {
		dRmin = dR;
		index = z;
	      }
	    }
	    
	    hlt_ph_offlineind[hlt_ph_n] = index;
	    hlt_ph_candpath[hlt_ph_n]= 0;
	    set(hlt_ph_candpath[hlt_ph_n], n);
	    hlt_ph_n++;
	  } else {
	    set(hlt_ph_candpath[isStored], n);
	  }
	}
      }
    }
    
    // MUON TRIGGER
    for(unsigned int n=0; n<theMuHLTLabels.size(); n++) {
      std::vector<edm::Ref<reco::RecoChargedCandidateCollection> > recoecalcands;
      trigger::size_type elIndex = triggerObjWithRef->filterIndex(theMuHLTLabels[n]);
      
      if (!(elIndex >= triggerObjWithRef->size())) {
	triggerObjWithRef->getObjects(elIndex, trigger::TriggerMuon, recoecalcands);
	
	for (unsigned int i=0; i<recoecalcands.size(); i++) {
	  if (hlt_el_n >= MAX_HLT) {
	    std::cout << "GlobeHLT: WARNING TOO MANY MUON HLT CANDIDATE: (allowed " << MAX_HLT << ")" << std::endl;
	    break;
	  }
	  
	  int isStored = -1;
	  TLorentzVector lv(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->energy());
	  for(int j=0; j<hlt_mu_n; j++) {
	    TLorentzVector* lv_temp =  (TLorentzVector*)hlt_mu_p4->At(j);
	    if (lv == *lv_temp) {
	      isStored = j;
	      break;
	    }
	  }
	  
	  if (isStored == -1) {
	    new ((*hlt_mu_p4)[hlt_mu_n]) TLorentzVector();
	    ((TLorentzVector *)hlt_mu_p4->At(hlt_mu_n))->SetXYZT(recoecalcands[i]->px(), recoecalcands[i]->py(), recoecalcands[i]->pz(), recoecalcands[i]->p());
	    
	    float pt = recoecalcands[i]->pt();
	    float eta = recoecalcands[i]->eta();
	    float phi = recoecalcands[i]->phi();
	    TVector3 v1;
	    v1.SetPtEtaPhi(pt, eta, phi);
	    
	    float dRmin = 1.0;
	    int index = -1;
	    for(unsigned int z=0; z<muH->size(); z++) {
	      reco::MuonRef e(muH, z);
	      // apply the cuts
	      if(gCUT->cut(*e))
		continue;
	      // passed cuts
	      
	      TVector3 v2;
	      v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	      float dR = v1.DeltaR(v2);
	      if (dR < dRmin) {
		dRmin = dR;
		index = z;
	      }
	    }
	    
	    hlt_mu_offlineind[hlt_mu_n] = index;
	    hlt_mu_candpath[hlt_mu_n]= 0;
	    set(hlt_mu_candpath[hlt_mu_n], n);
	    hlt_mu_n++;
	  } else {
	    set(hlt_mu_candpath[isStored], n);
	  }
	}
      }
    }
    ////////
  } else {
    ////////
    // HLT1
    hlt1_bit_1 = 0;
    hlt1_bit_2 = 0;
    hlt1_bit_3 = 0;
    hlt1_bit_4 = 0;
    bool changed = false;
    configProvider.init(iEvent.getRun(),iSetup,hlt1Tag_.process(),changed);

    edm::Handle<edm::TriggerResults> h_triggerResults_HLT1;
    iEvent.getByLabel(hlt1Tag_, h_triggerResults_HLT1);
    if (h_triggerResults_HLT1.isValid()) {
      hlt_path_names_HLT1_1->clear();
      hlt_path_names_HLT1_2->clear();
      hlt_path_names_HLT1_3->clear();
      hlt_path_names_HLT1_4->clear();
      if(debug_level > 9) std::cout << "Fill names HLT1" << std::endl;
      for (size_t i = 0; i < configProvider.size(); ++i) {
	int j=(int)(i/32);
	if(j<1) {     
	  hlt_path_names_HLT1_1->push_back(configProvider.triggerName(i));
	} else if(j<2){
	  hlt_path_names_HLT1_2->push_back(configProvider.triggerName(i));
	} else if(j<3){
	    hlt_path_names_HLT1_3->push_back(configProvider.triggerName(i));
	} else if(j<4){
	  hlt_path_names_HLT1_4->push_back(configProvider.triggerName(i));
	}
      }
      // Trigger Results
      if(debug_level > 99) std::cout << "### Trigger Results 1 :" << hlt1Tag_.process() << std::endl;
      for (size_t i = 0; i < configProvider.size(); ++i) {
	if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT1->accept(i) ? "passed" : "failed") << std::endl;
	int j=(int)(i/32);
	if(j<1) {     
	  if(h_triggerResults_HLT1->accept(i))
	    set(hlt1_bit_1,i-32*j);
	} else if(j<2){
	  if(h_triggerResults_HLT1->accept(i))
	    set(hlt1_bit_2,i-32*j);
	} else if(j<3){
	  if(h_triggerResults_HLT1->accept(i))
	    set(hlt1_bit_3,i-32*j);
	} else if(j<4){
	  if(h_triggerResults_HLT1->accept(i))
	    set(hlt1_bit_4,i-32*j);
	}
	//
	if(debug_level > 999) {
	  std::bitset<32> binary_hlt1(hlt1_bit_1);
	  std::bitset<32> binary_hlt2(hlt1_bit_2);
	  std::bitset<32> binary_hlt3(hlt1_bit_3);
	  std::bitset<32> binary_hlt4(hlt1_bit_4);
	  std::cout << "HLT Path 1: " << binary_hlt1  << " ? "<< check(hlt1_bit_1,i-32*j) << std::endl;
	  std::cout << "HLT Path 2: " << binary_hlt2  << " ? "<< check(hlt1_bit_2,i-32*j) << std::endl;
	  std::cout << "HLT Path 3: " << binary_hlt3  << " ? "<< check(hlt1_bit_3,i-32*j) << std::endl;
	  std::cout << "HLT Path 4: " << binary_hlt4  << " ? "<< check(hlt1_bit_4,i-32*j) << std::endl;
	}
	//
      }
      if(debug_level > 99) std::cout << "\t Final result = " << hlt1_bit_1 << " " << hlt1_bit_2 << " " << hlt1_bit_3 << " " << hlt1_bit_4 << std::endl;
    } else {
      if(debug_level > 9) std::cout << "TriggerResults not valid " << hlt1Tag_ << std::endl;
    }
    
    // HLT2
    if(secondaryTriggerON){
      hlt2_bit_1 = 0;
      hlt2_bit_2 = 0;
      hlt2_bit_3 = 0;
      hlt2_bit_4 = 0;
      bool changed = false;
      configProvider.init(iEvent.getRun(),iSetup,hlt2Tag_.process(),changed);
      edm::Handle<edm::TriggerResults> h_triggerResults_HLT2;
      iEvent.getByLabel(hlt2Tag_, h_triggerResults_HLT2);
      if (h_triggerResults_HLT2.isValid()) {
        hlt_path_names_HLT2_1->clear();
        hlt_path_names_HLT2_2->clear();
        hlt_path_names_HLT2_3->clear();
        hlt_path_names_HLT2_4->clear();
        if(debug_level > 9) std::cout << "Fill names HLT2" << std::endl;
        for (size_t i = 0; i < configProvider.size(); ++i) {
          int j=(int)(i/32);
          if(j<1) {     
            hlt_path_names_HLT2_1->push_back(configProvider.triggerName(i));
          } else if(j<2){
            hlt_path_names_HLT2_2->push_back(configProvider.triggerName(i));
          } else if(j<3){
            hlt_path_names_HLT2_3->push_back(configProvider.triggerName(i));
          } else if(j<4){
            hlt_path_names_HLT2_4->push_back(configProvider.triggerName(i));
          }
        }
        // Trigger Results
        if(debug_level > 99) std::cout << "### Trigger Results 2: " << hlt2Tag_.process() << std::endl;
        for (size_t i = 0; i < configProvider.size(); ++i) {
          if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT2->accept(i) ? "passed" : "failed") << std::endl;
          int j=(int)(i/32);
          if(j<1) {     
            if(h_triggerResults_HLT2->accept(i))
              set(hlt2_bit_1,i-32*j);
          } else if(j<2){
            if(h_triggerResults_HLT2->accept(i))
              set(hlt2_bit_2,i-32*j);
          } else if(j<3){
            if(h_triggerResults_HLT2->accept(i))
              set(hlt2_bit_3,i-32*j);
          } else if(j<4){
            if(h_triggerResults_HLT2->accept(i))
              set(hlt2_bit_4,i-32*j);
          }
          //
          if(debug_level > 999) {
            std::bitset<32> binary_hlt1(hlt2_bit_1);
            std::bitset<32> binary_hlt2(hlt2_bit_2);
            std::bitset<32> binary_hlt3(hlt2_bit_3);
            std::bitset<32> binary_hlt4(hlt2_bit_4);
            std::cout << "HLT Path 1: " << binary_hlt1  << " ? "<< check(hlt2_bit_1,i-32*j) << std::endl;
            std::cout << "HLT Path 2: " << binary_hlt2  << " ? "<< check(hlt2_bit_2,i-32*j) << std::endl;
            std::cout << "HLT Path 3: " << binary_hlt3  << " ? "<< check(hlt2_bit_3,i-32*j) << std::endl;
            std::cout << "HLT Path 4: " << binary_hlt4  << " ? "<< check(hlt2_bit_4,i-32*j) << std::endl;
          }
          //
        }
        if(debug_level > 99) std::cout << "\t Final result = " << hlt2_bit_1 << " " << hlt2_bit_2 << " " << hlt2_bit_3 << " " << hlt2_bit_4 << std::endl;
      } else {
        if(debug_level > 9) std::cout << "TriggerResults not valid " << hlt2Tag_ << std::endl;
      }
    }
    
    if(!triggerObj.isValid()) 
      throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );
    // Do it only the first event (when theHLTLabels is still empty)
    if(debug_level > 9) std::cout << "Trigger Paths" << std::endl;
    hlt_label_names_1->clear();
    hlt_label_names_2->clear();
    hlt_label_names_3->clear();
    hlt_label_names_4->clear();
    theHLTLabels.clear();
    for(int i=0; i<triggerObj->sizeFilters(); ++i) {
      if(debug_level > 9) std::cout << triggerObj->filterTag(i) << std::endl;
      theHLTLabels.push_back( triggerObj->filterTag(i));
      int j=(int)(i/32);
      if(j<1) {
	hlt_label_names_1->push_back(triggerObj->filterTag(i).label());
      } else if(j<2) {
	hlt_label_names_2->push_back(triggerObj->filterTag(i).label());
      } else if(j<3) {
	hlt_label_names_3->push_back(triggerObj->filterTag(i).label());
      } else if(j<4) {
	hlt_label_names_4->push_back(triggerObj->filterTag(i).label());
      }
    }
    
    hlt_n = 0;
    hlt_p4->Clear();
    
    trigger::TriggerObjectCollection triggerObjs = triggerObj->getObjects();
    if(debug_level > 99) std::cout << "Trigger Objects found " << triggerObjs.size() << std::endl;
    
    for (unsigned int iCand=0; iCand<triggerObjs.size(); ++iCand ) {
      if(debug_level > 99) std::cout << iCand << "=" << hlt_n << std::endl;
      if (hlt_n >= MAX_HLT) {
	std::cout << "GlobeHLT: WARNING TOO MANY ELECTRON HLT CANDIDATE: (allowed " << MAX_HLT << ")" << std::endl;
	break;
      }
      
      trigger::TriggerObject object = triggerObjs[iCand];
      hlt_candpath_1[hlt_n] = 0;
      hlt_candpath_2[hlt_n] = 0;
      hlt_candpath_3[hlt_n] = 0;
      hlt_candpath_4[hlt_n] = 0;
      
      for(unsigned int n=0; n<theHLTLabels.size(); n++) {
	//std::vector<edm::Ref<reco::ElectronCollection> > recoecalcands;
	trigger::size_type elIndex = triggerObj->filterIndex(theHLTLabels[n]);
	// std::cout << "INDEX: " << elIndex << std::endl;
	
	// Check HLT
	bool firedHLT = false;
	if (!(elIndex >= triggerObj->sizeFilters())) {
	  const trigger::Keys & k = triggerObj->filterKeys(elIndex);
	  for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
	    if(*ki == iCand)
	      firedHLT = true;
	  }
	}
	
	
	if(firedHLT) {
	  if(debug_level > 99) std::cout << n << "\t" << theHLTLabels[n].label() << "\t fired" << std::endl;
	  int j=(int)(n/32);
	  if(j<1) {
	    set(hlt_candpath_1[hlt_n], n-32*j);
	  } else if(j<2) {
	    set(hlt_candpath_2[hlt_n], n-32*j);
	  } else if(j<3) {
	    set(hlt_candpath_3[hlt_n], n-32*j);
	  } else if(j<4) {
	    set(hlt_candpath_4[hlt_n], n-32*j);
	  }
	}
	
      }
      
      // Skip if no triggers were fired
      if(hlt_candpath_1[hlt_n]!=0 || hlt_candpath_2[hlt_n]!=0 || hlt_candpath_3[hlt_n]!=0 || hlt_candpath_4[hlt_n]!=0) {
	
	// Set variables
	hlt_id[hlt_n] = object.id();
	
	// Set HLT candidate p4
	TLorentzVector lv(object.px(), object.py(), object.pz(), 0);
	new ((*hlt_p4)[hlt_n]) TLorentzVector();
	((TLorentzVector *)hlt_p4->At(hlt_n))->SetXYZT(object.px(), object.py(), object.pz(), object.energy());
	
	// OFFLINE matchings           
	float pt = object.pt(); 
	float eta = object.eta(); 
	float phi = object.phi(); 
	TVector3 v1;
	v1.SetPtEtaPhi(pt, eta, phi);
	// ELECTRON Matching
	float dRmin = 1.0;
	int index = -1;
	for(unsigned int z=0; z<elH->size(); z++) {
	  reco::GsfElectronRef e(elH, z);
	  // apply the cuts
	  if(gCUT->cut(*e))
	    continue;
	  // passed cuts
	  TVector3 v2;
	  v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	  float dR = v1.DeltaR(v2); 
	  if (dR < dRmin) {
	    dRmin = dR;
	    index = z;
	  }
	}
	hlt_el_offlineind[hlt_n] = index;
	
	// PHOTON matching
	dRmin = 1.0;
	index = -1;
	for(unsigned int z=0; z<phH->size(); z++) {
	  reco::PhotonRef e(phH, z);
	  // apply the cuts
	  if(gCUT->cut(*e))
	    continue;
	  // passed cuts
	  TVector3 v2;
	  v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	  float dR = v1.DeltaR(v2); 
	  if (dR < dRmin) {
	    dRmin = dR;
	    index = z;
	  }
	}
	hlt_ph_offlineind[hlt_n] = index;
	
	// MUON matching
	dRmin = 1.0;
	index = -1;
	for(unsigned int z=0; z<muH->size(); z++) {
	  reco::MuonRef e(muH, z);
	  // apply the cuts
	  if(gCUT->cut(*e))
	    continue;
	  // passed cuts
	  TVector3 v2;
	  v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	  float dR = v1.DeltaR(v2); 
	  if (dR < dRmin) {
	    dRmin = dR;
	    index = z;
	  }
	}
	hlt_mu_offlineind[hlt_n] = index; 
	
	// JET matching
	dRmin = 1.0;
	index = -1;
	for(unsigned int z=0; z<jetH->size(); z++) {
	  reco::CaloJetRef e(jetH, z);
	  // apply the cuts
	  if(gCUT->cut(*e))
	    continue;
	  // passed cuts
	  TVector3 v2;
	  v2.SetPtEtaPhi(e->pt(), e->eta(), e->phi());
	  float dR = v1.DeltaR(v2); 
	  if (dR < dRmin) {
	    dRmin = dR;
	    index = z;
	  }
	}
	hlt_jet_offlineind[hlt_n] = index; 
	
	//
	if(debug_level > 99) {
	  std::bitset<32> binary_hlt1(hlt_candpath_1[hlt_n]);
	  std::bitset<32> binary_hlt2(hlt_candpath_2[hlt_n]);
	  std::bitset<32> binary_hlt3(hlt_candpath_3[hlt_n]);
	  std::bitset<32> binary_hlt4(hlt_candpath_4[hlt_n]);
	  std::cout << iCand << "=" << hlt_n
		    <<"\t pT=" << object.pt() << "\t eta=" << object.eta() << "\t particle " << object.id() << std::endl
		    << "\t candpath1=" << binary_hlt1 << "\t candpath2=" << binary_hlt2
		    << "\t candpath3=" << binary_hlt3 << "\t candpath4=" << binary_hlt4
		    << "\t Matching"
		    << "\t electron " << hlt_el_offlineind[ hlt_n]
		    << "\t photon "   << hlt_ph_offlineind[ hlt_n]
		    << "\t muon "     << hlt_mu_offlineind[ hlt_n]
		    << "\t jet "      << hlt_jet_offlineind[hlt_n]
		    << std::endl;
	}
	//
	hlt_n++;          
      } // Store candidate which fired at least 1 HLT
      
    } // TriggerCandidate's Loop
    
    if(debug_level > 99) std::cout << "Trigger Objects stored " << hlt_n << std::endl;
  } // full HLT or not
  
  return true;
}

// USEFUL LINKS
// https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_8E29
// http://cms-project-confdb-hltdev.web.cern.ch/cms-project-confdb-hltdev/browser/convert2Html.jsp?dbName=HLTDEV&configName=/CMSSW/CMSSW_3_1_2/8E29/V1
// https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_1E31
// http://cms-project-confdb-hltdev.web.cern.ch/cms-project-confdb-hltdev/browser/convert2Html.jsp?dbName=HLTDEV&configName=/CMSSW/CMSSW_3_1_2/1E31/V1
// http://cms-service-sdtweb.web.cern.ch/cms-service-sdtweb/doxygen/CMSSW_3_2_5/doc/html/da/d01/TriggerTypeDefs_8h-source.html

