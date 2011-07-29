#ifndef GLOBEJETS_H
#define GLOBEJETS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeJets {
 public:
  
  GlobeJets(const edm::ParameterSet&, const char*);
  virtual ~GlobeJets() {};

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
 
  // variables
  Int_t jet_n;

  Float_t jet_emfrac[MAX_JETS];
  Float_t jet_hadfrac[MAX_JETS];

  Float_t jet_pull_dy[MAX_JETS];
  Float_t jet_pull_dphi[MAX_JETS];
  Float_t jet_erescale[MAX_JETS];

  Int_t jet_ntk[MAX_JETS];
  Int_t jet_ncalotw[MAX_JETS];
  std::vector<std::vector<unsigned short> >* jet_calotwind;
  std::vector<std::vector<unsigned short> >* jet_tkind;
  
  TClonesArray *jet_p4;

 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag jetColl, calotowerColl, trackColl, jetTkAssColl;
  std::string pfak5corrdata, pfak5corrmc, pfak5corr;
  edm::InputTag bcBColl, bcEColl, tkColl, pfJetColl;
  int debug_level;

  bool doEgammaSummer09Skim;
};

#endif
