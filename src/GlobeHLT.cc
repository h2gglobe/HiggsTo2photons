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
}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  //edm::Handle<trigger::TriggerEventWithRefs> triggerObjWithRef;
  
  iEvent.getByLabel(inputTag_, triggerObj);
  
  // HLT
  hlt_bit->clear();

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
      //if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT1->accept(i) ? "passed" : "failed") << std::endl;
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
  
  return true;
}
