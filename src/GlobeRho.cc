#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeRho.h"

GlobeRho::GlobeRho(const edm::ParameterSet& iConfig, const char* n):nome(n) {
  
  char a[100];
  sprintf(a, "rhoCollection_%s", nome);
  rhoCollection =  iConfig.getParameter<edm::InputTag>(a);
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeRho::defineBranch(TTree* tree) {
  
  char a1[100], a2[100];
  sprintf(a1, "rho_%s", nome);
  sprintf(a2, "rho_%s/F", nome);
  tree->Branch(a1, &rho, a2);
}

bool GlobeRho::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoCollection, rhoHandle);
  rho = *(rhoHandle.product());
 
  return true;
}
