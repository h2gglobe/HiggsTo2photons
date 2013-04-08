#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeRho.h"

GlobeRho::GlobeRho(const edm::ParameterSet& iConfig) {
  
  rhoCollection =  iConfig.getParameter<edm::InputTag>("rhoCorrection");
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeRho::defineBranch(TTree* tree) {
  
  tree->Branch("rho", &rho,"rho/F");
}

bool GlobeRho::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  if(strcmp(rhoCollection.encode().c_str(), "") != 0) {
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(rhoCollection, rhoHandle);
    rho = *(rhoHandle.product());
  } else {
    rho = 99999.;
  }

  return true;
}
