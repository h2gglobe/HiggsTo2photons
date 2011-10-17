#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenMET.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

GlobeGenMET::GlobeGenMET(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  METColl =  iConfig.getParameter<edm::InputTag>("GenTrueMETColl");

  if(strcmp(nome, "calo") == 0) 
    METColl =  iConfig.getParameter<edm::InputTag>("GenCaloMETColl");
  
  if(strcmp(nome, "nopt") == 0)
    METColl =  iConfig.getParameter<edm::InputTag>("GenNoptMETColl");
  
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeGenMET::defineBranch(TTree* tree) {

  char a1[100], a2[100];
  sprintf(a1, "met_%s_met", nome);
  sprintf(a2, "met_%s_met/F", nome);
  tree->Branch(a1, &met_met, a2);

  sprintf(a1, "met_%s_phi", nome);
  sprintf(a2, "met_%s_phi/F", nome);
  tree->Branch(a1, &met_phi, a2);
}

bool GlobeGenMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  if(debug_level > 99 ) 
    std::cout << "GlobeGenMET: Start " << std::endl;
    
  edm::Handle<reco::GenMETCollection> metH;
  iEvent.getByLabel(METColl, metH);
    
  met_met = metH->begin()->et();
  met_phi = metH->begin()->phi();
        
  if(debug_level > 99 ) 
    std::cout << "GlobeGenMET: End " << std::endl;
  
  return true;
}

