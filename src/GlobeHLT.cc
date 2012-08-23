#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHLT.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

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
  inputTag_ = psetHLT.getParameter<edm::InputTag>("TriggerResultsTag");
  hltTag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");

  hlt_path_names_HLT = new std::vector<std::string>;
  hlt_bit = new std::vector<unsigned short>;
  hlt_candpath = new std::vector<std::vector<unsigned short> >; 
  hlt_candpath2 = new std::vector<std::vector<unsigned short> >; 

  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeHLT::defineBranch(TTree* tree) {

  hlt_p4  = new TClonesArray("TLorentzVector", MAX_HLT);

  // Event Trigger
  tree->Branch("hlt_bit", "std::vector<unsigned short>", &hlt_bit);
  tree->Branch("hlt_path_names_HLT", "std::vector<std::string>", &hlt_path_names_HLT);

  // Trigger Candidates
  tree->Branch("hlt_n", &hlt_n, "hlt_n/I");
  tree->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0);
  tree->Branch("hlt_candpath", "std::vector<std::vector<unsigned short> >", &hlt_candpath);
  tree->Branch("hlt_candpath2", "std::vector<std::vector<unsigned short> >", &hlt_candpath2);

  //filter passes
  tree->Branch("filter_names_HLT1", "std::vector<std::string>", &filter_names_HLT1);
  tree->Branch("filter_pass", "std::vector<unsigned int>", &filter_pass);

  //Mass filter decisions
  //RELATED to making sure OR is same as four separate paths
  /*
  tree->Branch("pass_Mass60_isoiso",&pass_Mass60_isoiso, "pass_Mass60_isoiso/I");
  tree->Branch("pass_Mass60_R9R9",&pass_Mass60_R9R9, "pass_Mass60_R9R9/I");
  tree->Branch("pass_Mass60_mix",&pass_Mass60_mix, "pass_Mass60_mix/I");
  tree->Branch("pass_Mass70_isoiso",&pass_Mass70_isoiso, "pass_Mass70_isoiso/I");
  tree->Branch("pass_Mass70_R9R9",&pass_Mass70_R9R9, "pass_Mass70_R9R9/I");
  tree->Branch("pass_Mass70_mix",&pass_Mass70_mix, "pass_Mass70_mix/I");
*/

  //Ele32_SC17 trigger objects
  tree->Branch("trg_SC_ele_n", &ElectronRefs0_n,"ElectronRefs0_n/I");
  tree->Branch("trg_SC_ele_eta", &ElectronRefs0_eta,"ElectronRefs0_eta[ElectronRefs0_n]/F");
  tree->Branch("trg_SC_ele_et", &ElectronRefs0_et,"ElectronRefs0_et[ElectronRefs0_n]/F");
  tree->Branch("trg_SC_ele_phi", &ElectronRefs0_phi,"ElectronRefs0_phi[ElectronRefs0_n]/F");
  tree->Branch("trg_ele_n", &ElectronRefs1_n,"ElectronRefs1_n/I");
  tree->Branch("trg_ele_eta", &ElectronRefs1_eta,"ElectronRefs1_eta[ElectronRefs1_n]/F");
  tree->Branch("trg_ele_et", &ElectronRefs1_et,"ElectronRefs1_et[ElectronRefs1_n]/F");
  tree->Branch("trg_ele_phi", &ElectronRefs1_phi,"ElectronRefs1_phi[ElectronRefs1_n]/F");

  //26_18 OR trigger objects
  tree->Branch("PhotonRefs0_n", &PhotonRefs0_n,"PhotonRefs0_n/I");
  tree->Branch("PhotonRefs0_eta", &PhotonRefs0_eta,"PhotonRefs0_eta[PhotonRefs0_n]/F");
  tree->Branch("PhotonRefs0_et", &PhotonRefs0_et,"PhotonRefs0_et[PhotonRefs0_n]/F");
  tree->Branch("PhotonRefs0_phi", &PhotonRefs0_phi,"PhotonRefs0_phi[PhotonRefs0_n]/F");
  tree->Branch("PhotonRefs1_n", &PhotonRefs1_n,"PhotonRefs1_n/I");
  tree->Branch("PhotonRefs1_eta", &PhotonRefs1_eta,"PhotonRefs1_eta[PhotonRefs1_n]/F");
  tree->Branch("PhotonRefs1_et", &PhotonRefs1_et,"PhotonRefs1_et[PhotonRefs1_n]/F");
  tree->Branch("PhotonRefs1_phi", &PhotonRefs1_phi,"PhotonRefs1_phi[PhotonRefs1_n]/F");
  tree->Branch("PhotonRefs3_n", &PhotonRefs3_n,"PhotonRefs3_n/I");
  tree->Branch("PhotonRefs3_eta", &PhotonRefs3_eta,"PhotonRefs3_eta[PhotonRefs3_n]/F");
  tree->Branch("PhotonRefs3_et", &PhotonRefs3_et,"PhotonRefs3_et[PhotonRefs3_n]/F");
  tree->Branch("PhotonRefs3_phi", &PhotonRefs3_phi,"PhotonRefs3_phi[PhotonRefs3_n]/F");
  tree->Branch("PhotonRefs4_n", &PhotonRefs4_n,"PhotonRefs4_n/I");
  tree->Branch("PhotonRefs4_eta", &PhotonRefs4_eta,"PhotonRefs4_eta[PhotonRefs4_n]/F");
  tree->Branch("PhotonRefs4_et", &PhotonRefs4_et,"PhotonRefs4_et[PhotonRefs4_n]/F");
  tree->Branch("PhotonRefs4_phi", &PhotonRefs4_phi,"PhotonRefs4_phi[PhotonRefs4_n]/F");

  //36_22 OR trigger objects
  tree->Branch("PhotonRefs5_n", &PhotonRefs5_n,"PhotonRefs5_n/I");
  tree->Branch("PhotonRefs5_eta", &PhotonRefs5_eta,"PhotonRefs5_eta[PhotonRefs5_n]/F");
  tree->Branch("PhotonRefs5_et", &PhotonRefs5_et,"PhotonRefs5_et[PhotonRefs5_n]/F");
  tree->Branch("PhotonRefs5_phi", &PhotonRefs5_phi,"PhotonRefs5_phi[PhotonRefs5_n]/F");
  tree->Branch("PhotonRefs6_n", &PhotonRefs6_n,"PhotonRefs6_n/I");
  tree->Branch("PhotonRefs6_eta", &PhotonRefs6_eta,"PhotonRefs6_eta[PhotonRefs6_n]/F");
  tree->Branch("PhotonRefs6_et", &PhotonRefs6_et,"PhotonRefs6_et[PhotonRefs6_n]/F");
  tree->Branch("PhotonRefs6_phi", &PhotonRefs6_phi,"PhotonRefs6_phi[PhotonRefs6_n]/F");
  tree->Branch("PhotonRefs8_n", &PhotonRefs8_n,"PhotonRefs8_n/I");
  tree->Branch("PhotonRefs8_eta", &PhotonRefs8_eta,"PhotonRefs8_eta[PhotonRefs8_n]/F");
  tree->Branch("PhotonRefs8_et", &PhotonRefs8_et,"PhotonRefs8_et[PhotonRefs8_n]/F");
  tree->Branch("PhotonRefs8_phi", &PhotonRefs8_phi,"PhotonRefs8_phi[PhotonRefs8_n]/F");
  tree->Branch("PhotonRefs9_n", &PhotonRefs9_n,"PhotonRefs9_n/I");
  tree->Branch("PhotonRefs9_eta", &PhotonRefs9_eta,"PhotonRefs9_eta[PhotonRefs9_n]/F");
  tree->Branch("PhotonRefs9_et", &PhotonRefs9_et,"PhotonRefs9_et[PhotonRefs9_n]/F");
  tree->Branch("PhotonRefs9_phi", &PhotonRefs9_phi,"PhotonRefs9_phi[PhotonRefs9_n]/F");

  filter_pass = new std::vector<unsigned int>;
  filter_names_HLT1 = new std::vector<std::string>;
}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<trigger::TriggerEvent> triggerObj;

  iEvent.getByLabel(inputTag_, triggerObj);

  // HLT
  hlt_bit->clear();
  filter_pass->clear();
  filter_names_HLT1->clear();

  // FIXME possibly move to beginRun
  bool changed = false;
  configProvider.init(iEvent.getRun(), iSetup, hltTag_.process(), changed);
  edm::Handle<edm::TriggerResults> h_triggerResults_HLT;
  iEvent.getByLabel(hltTag_, h_triggerResults_HLT);

  if (h_triggerResults_HLT.isValid()) {
	hlt_path_names_HLT->clear();

	if(debug_level > 9) 
	  std::cout << "Fill names HLT" << std::endl;

	for (size_t i = 0; i < configProvider.size(); ++i)
	  hlt_path_names_HLT->push_back(configProvider.triggerName(i));

	// Trigger Results
	for (size_t i = 0; i < configProvider.size(); ++i) {
	  if(h_triggerResults_HLT->accept(i))
		hlt_bit->push_back((unsigned short)(i));
	}
  }

  if(!triggerObj.isValid()) 
	throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );

  hlt_n = 0;
  hlt_p4->Clear();
  hlt_candpath->clear();
  hlt_candpath2->clear();

  trigger::TriggerObjectCollection triggerObjects = triggerObj->getObjects();

  for (unsigned int iCand=0; iCand<triggerObjects.size(); ++iCand ) {

	if(debug_level > 99) 
	  std::cout << iCand << "=" << hlt_n << std::endl;
	if (hlt_n >= MAX_HLT) {
	  std::cout << "GlobeHLT: WARNING TOO MANY HLT CANDIDATES:  " << hlt_n << " found (allowed " << MAX_HLT << ")" << std::endl;
	  break;
	}

	std::vector<unsigned short> temp, temp2;
	trigger::TriggerObject object = triggerObjects[iCand];

	for (size_t i=0; i<configProvider.triggerNames().size(); i++) {
	  std::vector<std::string> labels = configProvider.saveTagsModules(configProvider.triggerNames()[i]); 
	  if (labels.size() == 0)
		continue;
	  edm::InputTag label(labels[labels.size()-1], "", hltTag_.process());
	  size_t filterIndex = triggerObj->filterIndex(label);
	  if (filterIndex < triggerObj->sizeFilters()) {
		const trigger::Keys &keys = triggerObj->filterKeys(filterIndex);
		for (size_t j = 0; j < keys.size(); j++) {
		  if (keys[j] == iCand)
			temp.push_back((unsigned short)i);
		}
	  } 

	  if (labels.size() > 2) {
		edm::InputTag label(labels[labels.size()-2], "", hltTag_.process());
		size_t filterIndex = triggerObj->filterIndex(label);

		if (filterIndex < triggerObj->sizeFilters()) {
		  const trigger::Keys &keys = triggerObj->filterKeys(filterIndex);
		  for (size_t j = 0; j < keys.size(); j++) {
			if (iCand == keys[j]) 
			  temp2.push_back((unsigned short)i);
		  }
		}
	  } 

	  if (temp.size() > temp2.size())
		temp2.push_back(999);

	}

	if (temp.size() > 0) {
	  hlt_candpath2->push_back(temp2);
	  hlt_candpath->push_back(temp);

	  // Set HLT candidate p4
	  TLorentzVector lv(object.px(), object.py(), object.pz(), 0);
	  new ((*hlt_p4)[hlt_n]) TLorentzVector();
	  ((TLorentzVector *)hlt_p4->At(hlt_n))->SetXYZT(object.px(), object.py(), object.pz(), object.energy());

	  hlt_n++;          
	}
  } 

  //added Trigger Objects
  edm::InputTag trigResultsTag("TriggerResults","","HLT");
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); 
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);

  std::vector<std::string> temp_names;
  temp_names.clear();
  temp_names.push_back("hltEG26CaloId10Iso50HcalIsoLastFilter");
  temp_names.push_back("hltEG26R9Id85LastFilter");
  //temp_names.push_back("HLTEG26R9Id85ORCaloId10Iso50LegCombLastFilter");
  temp_names.push_back("hltEG18R9Id85LastFilterUnseeded");
  temp_names.push_back("hltEG18CaloId10Iso50TrackIsoLastFilterUnseeded");

  temp_names.push_back("hltEG36CaloId10Iso50HcalIsoLastFilter");
  temp_names.push_back("hltEG36R9Id85LastFilter");
  //temp_names.push_back("HLTEG36R9Id85ORCaloId10Iso50LegCombLastFilter");
  temp_names.push_back("hltEG22R9Id85LastFilterUnseeded");
  temp_names.push_back("hltEG22CaloId10Iso50TrackIsoLastFilterUnseeded");

  //temp_names.push_back("hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");
  temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter");
  temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter");
  filter_pass->clear();
  filter_names_HLT1->clear();
  std::vector<std::string>::iterator filter_it;

  std::vector<trigger::TriggerObject> PhotonRefs0;
  std::vector<trigger::TriggerObject> PhotonRefs1;
  std::vector<trigger::TriggerObject> PhotonRefs3;
  std::vector<trigger::TriggerObject> PhotonRefs4;

  std::vector<trigger::TriggerObject> PhotonRefs5;
  std::vector<trigger::TriggerObject> PhotonRefs6;
  std::vector<trigger::TriggerObject> PhotonRefs8;
  std::vector<trigger::TriggerObject> PhotonRefs9;


//  std::vector<trigger::TriggerObject> PhotonRefs10;
  std::vector<trigger::TriggerObject> ElectronRefs0;
  std::vector<trigger::TriggerObject> ElectronRefs1;
  //std::vector<trigger::TriggerObject> ElectronRefs00;

  if ( trigEvent.isValid() ){
	for (filter_it = temp_names.begin(); filter_it != temp_names.end(); ++filter_it){
	  filter_names_HLT1->push_back((std::string)(*filter_it));
	  const trigger::TriggerObjectCollection & triggerObjects = trigEvent -> getObjects();
	  trigger::size_type filter1_idx = trigEvent -> filterIndex (edm::InputTag(*filter_it,"","HLT") ) ;   
	  trigger::size_type n_filters    = trigEvent -> sizeFilters();
	  if ( filter1_idx < n_filters ) {
		const trigger::Keys & triggerKeys ( trigEvent -> filterKeys ( filter1_idx ) );
		const int nkeys = triggerKeys.size();
		filter_pass->push_back(nkeys);
		for (int ikey = 0; ikey < nkeys; ++ikey ) {
		  if (*filter_it == "hltEG26CaloId10Iso50HcalIsoLastFilter") PhotonRefs0.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG26R9Id85LastFilter") PhotonRefs1.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG18R9Id85LastFilterUnseeded") PhotonRefs3.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG18CaloId10Iso50TrackIsoLastFilterUnseeded") PhotonRefs4.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter") ElectronRefs0.push_back(triggerObjects[ triggerKeys [ikey] ]);
		  else if (*filter_it == "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ElectronRefs1.push_back(triggerObjects[ triggerKeys [ikey]]);
		  //else if (*filter_it == "hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded") PhotonRefs10.push_back(triggerObjects[ triggerKeys [ ikey ] ]);

		  else if (*filter_it == "hltEG36CaloId10Iso50HcalIsoLastFilter") PhotonRefs5.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG36R9Id85LastFilter") PhotonRefs6.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG22R9Id85LastFilterUnseeded") PhotonRefs8.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		  else if (*filter_it == "hltEG22CaloId10Iso50TrackIsoLastFilterUnseeded") PhotonRefs9.push_back(triggerObjects[ triggerKeys [ ikey ] ]);
		}
	  }
	  else 
		filter_pass->push_back(0);
	}
  }

  //final mass filter with L1seeded and SC
  ElectronRefs0_n= 0;
  for (unsigned int i=0; i<ElectronRefs0.size(); i++) {
	if (ElectronRefs0_n >= 8) break;
	trigger::TriggerObject pho=ElectronRefs0[i];
	ElectronRefs0_eta[i] = pho.eta();
	ElectronRefs0_phi[i] = pho.phi();
	ElectronRefs0_et[i] = pho.et();
	ElectronRefs0_n++;
  }

  //L1seeded
  ElectronRefs1_n= 0;
  for (unsigned int i=0; i<ElectronRefs1.size(); i++) {
	if (ElectronRefs1_n >= 8) break;
	trigger::TriggerObject pho=ElectronRefs1[i];
	ElectronRefs1_eta[i] = pho.eta();
	ElectronRefs1_phi[i] = pho.phi();
	ElectronRefs1_et[i] = pho.et();
	ElectronRefs1_n++;
  }

  PhotonRefs0_n= 0;
  for (unsigned int i=0; i<PhotonRefs0.size(); i++) {
	if (PhotonRefs0_n >= 8) break;
	trigger::TriggerObject pho=PhotonRefs0[i];
	PhotonRefs0_eta[i] = pho.eta();
	PhotonRefs0_phi[i] = pho.phi();
	PhotonRefs0_et[i] = pho.et();
	PhotonRefs0_n++;
  }

  PhotonRefs1_n= 0;
  for (unsigned int i=0; i<PhotonRefs1.size(); i++) {
	if (PhotonRefs1_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs1[i];
	PhotonRefs1_eta[i] = pho.eta();
	PhotonRefs1_phi[i] = pho.phi();
	PhotonRefs1_et[i] = pho.et();
	PhotonRefs1_n++;
  }
  PhotonRefs3_n= 0;
  for(unsigned int i=0; i<PhotonRefs3.size(); i++) {
	if (PhotonRefs3_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs3[i];
	PhotonRefs3_eta[i] = pho.eta();
	PhotonRefs3_phi[i] = pho.phi();
	PhotonRefs3_et[i] = pho.et();
	PhotonRefs3_n++;
  }
  PhotonRefs4_n= 0;
  for(unsigned int i=0; i<PhotonRefs4.size(); i++) {
	if (PhotonRefs4_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs4[i];
	PhotonRefs4_eta[i] = pho.eta();
	PhotonRefs4_phi[i] = pho.phi();
	PhotonRefs4_et[i] = pho.et();
	PhotonRefs4_n++;
  }


  PhotonRefs5_n= 0;
  for (unsigned int i=0; i<PhotonRefs5.size(); i++) {
	if (PhotonRefs5_n >= 8) break;
	trigger::TriggerObject pho=PhotonRefs5[i];
	PhotonRefs5_eta[i] = pho.eta();
	PhotonRefs5_phi[i] = pho.phi();
	PhotonRefs5_et[i] = pho.et();
	PhotonRefs5_n++;
  }

  PhotonRefs6_n= 0;
  for (unsigned int i=0; i<PhotonRefs6.size(); i++) {
	if (PhotonRefs6_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs6[i];
	PhotonRefs6_eta[i] = pho.eta();
	PhotonRefs6_phi[i] = pho.phi();
	PhotonRefs6_et[i] = pho.et();
	PhotonRefs6_n++;
  }
  PhotonRefs8_n= 0;
  for(unsigned int i=0; i<PhotonRefs8.size(); i++) {
	if (PhotonRefs8_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs8[i];
	PhotonRefs8_eta[i] = pho.eta();
	PhotonRefs8_phi[i] = pho.phi();
	PhotonRefs8_et[i] = pho.et();
	PhotonRefs8_n++;
  }
  PhotonRefs9_n= 0;
  for(unsigned int i=0; i<PhotonRefs9.size(); i++) {
	if (PhotonRefs9_n >= 8) break;
	trigger::TriggerObject pho = PhotonRefs9[i];
	PhotonRefs9_eta[i] = pho.eta();
	PhotonRefs9_phi[i] = pho.phi();
	PhotonRefs9_et[i] = pho.et();
	PhotonRefs9_n++;
  }


  //RELATED to making sure OR is same as four separate paths
  /*
  //logic for Mass isoiso
  pass_Mass60_isoiso = 0;
  Mass60_isoiso = 0;
  pass_Mass70_isoiso = 0;

  if (PhotonRefs4.size() >=2) {
	TLorentzVector e1;  
	TLorentzVector e2;  
	TLorentzVector meson;
	float mass = 0.;
	for (int i = 0; i < (int)PhotonRefs4.size(); i++){
	  for (int j = i+1; j < (int)PhotonRefs4.size(); j++){
		e1.SetPtEtaPhiM(PhotonRefs4.at(i).pt(),PhotonRefs4.at(i).eta(),PhotonRefs4.at(i).phi(),0.); 
		e2.SetPtEtaPhiM(PhotonRefs4.at(j).pt(),PhotonRefs4.at(j).eta(),PhotonRefs4.at(j).phi(),0.); 
		meson = e1 + e2; 
		mass = meson.M(); 
		if (mass > Mass60_isoiso) Mass60_isoiso = mass;
		if (mass>60)  pass_Mass60_isoiso = 1;
		if (mass>70)  pass_Mass70_isoiso = 1;

	  }
	}	
  }     

  //logic for Mass R9R9
  pass_Mass60_R9R9 = 0;
  pass_Mass70_R9R9 = 0;

  if (PhotonRefs3.size() >=2) {
	TLorentzVector e1;  
	TLorentzVector e2;  
	TLorentzVector meson;
	float mass = 0.;
	for (int i = 0; i < (int)PhotonRefs3.size(); i++){
	  for (int j = i+1; j < (int)PhotonRefs3.size(); j++){
		e1.SetPtEtaPhiM(PhotonRefs3.at(i).pt(),PhotonRefs3.at(i).eta(),PhotonRefs3.at(i).phi(),0.); 
		e2.SetPtEtaPhiM(PhotonRefs3.at(j).pt(),PhotonRefs3.at(j).eta(),PhotonRefs3.at(j).phi(),0.); 
		meson = e1 + e2; 
		mass = meson.M();  
		if (mass>60)  pass_Mass60_R9R9 = 1; 
		if (mass>70)  pass_Mass70_R9R9 = 1; 
	  }
	}
  }
  //logic for Mass mix
  pass_Mass60_mix = 0;
  pass_Mass70_mix = 0;
  if (PhotonRefs3.size() >=1 && PhotonRefs4.size() >=1) {
	for (int i = 0; i < (int)PhotonRefs3.size(); i++){
	  for (int j = 0; j < (int)PhotonRefs4.size(); j++){
		TLorentzVector e1;  
		TLorentzVector e2; 	
		TLorentzVector meson;
		float mass = 0.;
		e1.SetPtEtaPhiM(PhotonRefs3.at(i).pt(),PhotonRefs3.at(i).eta(),PhotonRefs3.at(i).phi(),0); 
		e2.SetPtEtaPhiM(PhotonRefs4.at(j).pt(),PhotonRefs4.at(j).eta(),PhotonRefs4.at(j).phi(),0); 
		meson = e1 + e2; 
		mass = meson.M();  
		if (mass>60)  pass_Mass60_mix = 1; 
		if (mass>70)  pass_Mass70_mix = 1; 
	  }
	}
  }
  */ 

  return true;
}
