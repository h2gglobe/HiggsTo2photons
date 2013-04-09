#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePdfWeights.h"

#include <iostream>

GlobePdfWeights::GlobePdfWeights(const edm::ParameterSet& iConfig) {
  
  pdfweightsCollList = iConfig.getParameter<std::vector <edm::InputTag> >("PdfWeightsCollList");
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobePdfWeights::defineBranch(TTree* tree) {

 
  tree->Branch("weight_n", &weight_n, "weight_n/I");
  tree->Branch("pdf_weights", pdf_weights, "pdf_weights[weight_n]/F");

}

bool GlobePdfWeights::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  
  
  weight_n=0;
  
  for (unsigned int iPDF=0; iPDF<pdfweightsCollList.size(); iPDF++){
    edm::Handle<std::vector<double> > weightHandle;
    iEvent.getByLabel(pdfweightsCollList[iPDF], weightHandle);
    std::vector<double> weights = (*weightHandle);
    for (unsigned int j=0; j<weights.size(); j++){
      pdf_weights[weight_n]=weights[j];
      weight_n++;
    }
  }
  return true;
}
