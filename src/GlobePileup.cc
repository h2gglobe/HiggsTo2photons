#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePileup.h"

GlobePileup::GlobePileup(const edm::ParameterSet& iConfig) {
  
  pileupCollection =  iConfig.getParameter<edm::InputTag>("pileupInfoCollection");
  debug_level = iConfig.getParameter<int>("Debug_Level"); 
  
  pu_zpos = new std::vector<float>; pu_zpos->clear();
  pu_sumpt_lowpt = new std::vector<float>; pu_sumpt_lowpt->clear();
  pu_sumpt_highpt = new std::vector<float>; pu_sumpt_highpt->clear();
  pu_ntrks_lowpt = new std::vector<int>; pu_ntrks_lowpt->clear();
  pu_ntrks_highpt = new std::vector<int>; pu_ntrks_highpt->clear();
  
  h1 = new TH1D("pileup", "pileup", 100, 0, 100);
  h2 = new TH1D("pileup_true", "pileup_true", 1000, 0, 100);
}

void GlobePileup::defineBranch(TTree* tree) {
  
  tree->Branch("pu_n", &pu_n, "pu_n/I");
  tree->Branch("pu_n_true", &pu_n_true, "pu_n_true/F");
  tree->Branch("pu_bunchcrossing", &pu_bunchcrossing, "pu_bunchcrossing/I");
  tree->Branch("pu_zpos", "std::vector<float>", &pu_zpos);
  tree->Branch("pu_sumpt_lowpt", "std::vector<float>", &pu_sumpt_lowpt);
  tree->Branch("pu_sumpt_highpt", "std::vector<float>", &pu_sumpt_highpt);
  tree->Branch("pu_ntrks_lowpt", "std::vector<int>", &pu_ntrks_lowpt);
  tree->Branch("pu_ntrks_highpt", "std::vector<int>", &pu_ntrks_highpt);
}

bool GlobePileup::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  // New PU code that should work in 42X
  edm::Handle<std::vector< PileupSummaryInfo> > PupInfo;
  iEvent.getByLabel(pileupCollection, PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  pu_n = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    int pu_bunchcrossing = PVI->getBunchCrossing();
    if(pu_bunchcrossing == 0) {
      pu_n = PVI->getPU_NumInteractions();
      pu_n_true = PVI->getTrueNumInteractions();
      h1->Fill(pu_n);
      h2->Fill(pu_n_true);
      *pu_zpos = PVI->getPU_zpositions();
      *pu_sumpt_lowpt = PVI->getPU_sumpT_lowpT();
      *pu_sumpt_highpt = PVI->getPU_sumpT_highpT();
      *pu_ntrks_lowpt = PVI->getPU_ntrks_lowpT();
      *pu_ntrks_highpt = PVI->getPU_ntrks_highpT();
      
      break;
    }
  }
  
  return true;
}
