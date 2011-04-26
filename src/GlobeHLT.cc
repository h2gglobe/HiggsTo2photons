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
  secondaryTriggerON  = psetHLT.getParameter<bool>("useSecondaryTrigger");
  
  hlt1Tag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");
  hlt_path_names_HLT1 = new std::vector<std::string>; hlt_path_names_HLT1->clear();

  if(secondaryTriggerON){
    hlt2Tag_  = psetHLT.getParameter<edm::InputTag>("SecondaryTriggerResultsTag");
    hlt_path_names_HLT2 = new std::vector<std::string>; hlt_path_names_HLT2->clear();
  }
  
  hlt1_bit = new std::vector<unsigned short>; hlt1_bit->clear();
  hlt2_bit = new std::vector<unsigned short>; hlt2_bit->clear();
  hlt_candpath = new std::vector<std::vector<unsigned short> >; hlt_candpath->clear();
      
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeHLT::defineBranch(TTree* tree) {
  
  hlt_p4  = new TClonesArray("TLorentzVector",   MAX_HLT);
  
  // Event Trigger
  tree->Branch("hlt1_bit", "std::vector<unsigned short>", &hlt1_bit);
  tree->Branch("hlt_path_names_HLT1", "std::vector<std::string>", &hlt_path_names_HLT1);

  if(secondaryTriggerON){
    tree->Branch("hlt2_bit", "std::vector<unsigned short>", &hlt2_bit);
    tree->Branch("hlt_path_names_HLT2", "std::vector<std::string>", &hlt_path_names_HLT2);
  }

  // Trigger Candidates
  tree->Branch("hlt_n", &hlt_n, "hlt_n/I");
  tree->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0);

  tree->Branch("hlt_candpath", "std::vector<std::vector<unsigned short> >", &hlt_candpath);
}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  edm::Handle<trigger::TriggerEventWithRefs> triggerObjWithRef;
  
  iEvent.getByLabel(inputTag_, triggerObj);
  
  // HLT1
  hlt1_bit->clear();

  bool changed = false;
  configProvider.init(iEvent.getRun(),iSetup,hlt1Tag_.process(),changed);
  edm::Handle<edm::TriggerResults> h_triggerResults_HLT1;
  iEvent.getByLabel(hlt1Tag_, h_triggerResults_HLT1);

  if (h_triggerResults_HLT1.isValid()) {
    hlt_path_names_HLT1->clear();

    if(debug_level > 9) 
      std::cout << "Fill names HLT1" << std::endl;
    for (size_t i = 0; i < configProvider.size(); ++i)
      hlt_path_names_HLT1->push_back(configProvider.triggerName(i));
    
    // Trigger Results
    if(debug_level > 99) 
      std::cout << "### Trigger Results 1 :" << hlt1Tag_.process() << std::endl;
    for (size_t i = 0; i < configProvider.size(); ++i) {
      if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT1->accept(i) ? "passed" : "failed") << std::endl;

      if(h_triggerResults_HLT1->accept(i))
        hlt1_bit->push_back((unsigned short)(i));
    }
  }

  // HLT2
  if(secondaryTriggerON){
    hlt2_bit->clear();
    
    bool changed = false;
    configProvider.init(iEvent.getRun(),iSetup,hlt2Tag_.process(),changed);
    
    edm::Handle<edm::TriggerResults> h_triggerResults_HLT2;
    iEvent.getByLabel(hlt2Tag_, h_triggerResults_HLT2);
    
    if (h_triggerResults_HLT2.isValid()) {
      hlt_path_names_HLT2->clear();
      
      if(debug_level > 9) 
        std::cout << "Fill names HLT2" << std::endl;
      for (size_t i = 0; i < configProvider.size(); ++i)
        hlt_path_names_HLT2->push_back(configProvider.triggerName(i));
      
      // Trigger Results
      if(debug_level > 99) 
        std::cout << "### Trigger Results 2: " << hlt2Tag_.process() << std::endl;
      
      for (size_t i = 0; i < configProvider.size(); ++i) {
        if(debug_level > 99) std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT2->accept(i) ? "passed" : "failed") << std::endl;
        
        if(h_triggerResults_HLT2->accept(i))
          hlt2_bit->push_back(i);
        
      }
    }
  }

  if(!triggerObj.isValid()) 
    throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );

  // This can be improved doing it only when necessary
  theHLTLabels.clear();
  for(int i=0; i<triggerObj->sizeFilters(); ++i) {
    if(debug_level > 9) 
      std::cout << triggerObj->filterTag(i) << std::endl;
    theHLTLabels.push_back( triggerObj->filterTag(i));
  }
 
  hlt_n = 0;
  hlt_p4->Clear();
  hlt_candpath->clear();
  
  trigger::TriggerObjectCollection triggerObjs = triggerObj->getObjects();
  if(debug_level > 99) 
    std::cout << "Trigger Objects found " << triggerObjs.size() << std::endl;
  
  for (unsigned int iCand=0; iCand<triggerObjs.size(); ++iCand ) {
    if(debug_level > 99) 
      std::cout << iCand << "=" << hlt_n << std::endl;
    if (hlt_n >= MAX_HLT) {
      std::cout << "GlobeHLT: WARNING TOO MANY HLT CANDIDATES:  " << hlt_n << " found (allowed " << MAX_HLT << ")" << std::endl;
      break;
    }
    
    std::vector<unsigned short> temp;
    trigger::TriggerObject object = triggerObjs[iCand];

    for(unsigned int n=0; n<theHLTLabels.size(); n++) {

      trigger::size_type index = triggerObj->filterIndex(theHLTLabels[n]);
      
      // Check HLT
      bool firedHLT = false;
      if (!(index >= triggerObj->sizeFilters())) {
        const trigger::Keys & k = triggerObj->filterKeys(index);
        for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
          if(*ki == iCand)
            firedHLT = true;
        }
      }
      
      
      if(firedHLT) {
        if(debug_level > 99) std::cout << n << "\t" << theHLTLabels[n].label() << "\t fired" << std::endl;
        
        for (unsigned int i=0; i<configProvider.size(); i++) {
          
          unsigned int nModules = configProvider.moduleLabels(i).size();
          unsigned int moduleIndex = configProvider.moduleIndex(i, theHLTLabels[n].label());
          //std::cout << nModules << " " << moduleIndex << std::endl;
          //for(unsigned int y=0; y<nModules; y++) {
          //  std::cout << y << " " << configProvider.moduleLabel(i, y) << " " << theHLTLabels[n].label() << " " << configProvider.triggerName(i) << std::endl;
          //  if (configProvider.moduleLabel(i, y) == theHLTLabels[n].label()) 
          //    std::cout << "MATCH !!!!"  << " " << y << " " << nModules << std::endl;
          // }
          //std::cout << "------------" << std::endl;
          if ((nModules - moduleIndex) == 2) {
            //std::cout <<  theHLTLabels[n].label() << " " << std::endl;;
            //std::cout << configProvider.moduleLabel(i, moduleIndex) << " " << nModules << " " << moduleIndex << std::endl;
            //std::cout <<  configProvider.triggerName(i) << std::endl;
            //std::cout << std::endl;
            temp.push_back(i);
          }
        }
      }
    }

    // Skip if no triggers were fired
    if(temp.size() != 0) {
      hlt_candpath->push_back(temp);
      
      // Set HLT candidate p4
      TLorentzVector lv(object.px(), object.py(), object.pz(), 0);
      new ((*hlt_p4)[hlt_n]) TLorentzVector();
      ((TLorentzVector *)hlt_p4->At(hlt_n))->SetXYZT(object.px(), object.py(), object.pz(), object.energy());
      
      hlt_n++;          
    } // Store candidate which fired at least 1 HLT
    
  } // TriggerCandidate's Loop
  
  if(debug_level > 99) std::cout << "Trigger Objects stored " << hlt_n << std::endl;
  
  return true;
}
