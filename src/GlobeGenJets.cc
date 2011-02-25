#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenJets.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <iostream>

GlobeGenJets::GlobeGenJets(const edm::ParameterSet& iConfig, const char* n = "algo1"): nome(n) {
  
  char a[100];
  sprintf (a,"GenJetColl_%s", nome);
  jetColl =  iConfig.getParameter<edm::InputTag>(a);
  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeGenJets::defineBranch(TTree* tree) {

  genjet_p4 = new TClonesArray("TLorentzVector", MAX_GENJETS);
  
  char a1[50], a2[50];
  
  sprintf(a1, "genjet_%s_n", nome);
  sprintf(a2, "genjet_%s_n/I", nome);
  tree->Branch(a1, &genjet_n, a2);
  
  sprintf(a1, "genjet_%s_p4", nome);
  tree->Branch(a1, "TClonesArray", &genjet_p4, 32000, 0);
  
  sprintf(a1, "genjet_%s_em", nome);
  sprintf(a2, "genjet_%s_em[genjet_%s_n]/F", nome, nome);
  tree->Branch(a1, &genjet_em, a2);

  sprintf(a1, "genjet_%s_had", nome);
  sprintf(a2, "genjet_%s_had[genjet_%s_n]/F", nome, nome);
  tree->Branch(a1, &genjet_had, a2);

  sprintf(a1, "genjet_%s_inv", nome);
  sprintf(a2, "genjet_%s_inv[genjet_%s_n]/F", nome, nome);
  tree->Branch(a1, &genjet_inv, a2);

  sprintf(a1, "genjet_%s_aux", nome);
  sprintf(a2, "genjet_%s_aux[genjet_%s_n]/F", nome, nome);
  tree->Branch(a1, &genjet_aux, a2);
}

bool GlobeGenJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // take collections
  edm::Handle<reco::GenJetCollection> jetH;
  iEvent.getByLabel(jetColl, jetH);

  genjet_p4->Clear();
  
  genjet_n = 0;
  if (debug_level > 9)
    std::cout << "GlobeGenJets: GenJet collection size: "<< jetH->size() << std::endl;
  
  // check if collection is present
  for(unsigned int i=0; i<jetH->size(); i++) {
    if (genjet_n >= MAX_GENJETS) {
      std::cout << "GlobeGenJets: WARNING TOO MANY GENJETS: " << jetH->size() << " (allowed " << MAX_GENJETS << ")" << std::endl;
      break;
    }
    reco::GenJetRef j(jetH, i);
	 // apply the cuts
	 if(gCUT->cut(*j))continue;
	 // passed cuts

    new ((*genjet_p4)[genjet_n]) TLorentzVector();
    ((TLorentzVector *)genjet_p4->At(genjet_n))->SetXYZT(j->px(), j->py(), j->pz(), j->energy()); 
    genjet_em[genjet_n] = j->emEnergy();
    genjet_had[genjet_n] = j->hadEnergy();
    genjet_inv[genjet_n] = j->invisibleEnergy();
    genjet_aux[genjet_n] = j->auxiliaryEnergy();
    genjet_n++;
    
  }
  
  //if (genjet_n == 0)
  //  return false;
  
  return true;
}
