#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEREDUCEDGEN_H
#define GLOBEREDUCEDGEN_H

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenerator.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <iostream>

class GlobeGenerator;
class GlobeGenParticles;
class GlobeAnalyzer;

class GlobeReducedGen {
 public:
  
  GlobeReducedGen(const edm::ParameterSet&);
  virtual ~GlobeReducedGen() {};
 
  void defineBranch(GlobeAnalyzer* ana);
  void fillRedGenList(GlobeGenerator * gen, GlobeLeptons * lep);
  int  LeptonsGenInfo(GlobeGenerator * gen, int j, int * gen_keep, int & nhistory, int * history, int & nrad_photons, int * rad_photons);
  int GlobeMatchWithGen(GlobeGenerator * gen, TLorentzVector* p4, Float_t & deltaR, int pdgid, int gencoll, float cutgenrecdrmatch);
  void fillRedGenList(GlobeGenParticles* gen, GlobeLeptons * lep);
  int  LeptonsGenInfo(GlobeGenParticles* gen, int j, int * gen_keep, int & nhistory, int * history, int & nrad_photons, int * rad_photons);
  int GlobeMatchWithGen(GlobeGenParticles* gen, TLorentzVector* p4, Float_t & deltaR, int pdgid, int gencoll, float cutgenrecdrmatch);

  // bool analyze(const edm::Event&, const edm::EventSetup&);
  //  int mother(const HepMC::GenEvent* g, HepMC::GenParticle *p);

  // variables
  Int_t lptgeninfo_n;
  Int_t lptgen_n;
  
  Int_t lptgeninfo_status[MAX_LPT_GENINFO];
  Int_t lptgeninfo_pdgid[MAX_LPT_GENINFO];
  Int_t lptgeninfo_mother[MAX_LPT_GENINFO];
  
  Int_t  lptgen_status[MAX_LPT_GENINFO];
  Int_t lptgen_pdgid[MAX_LPT_GENINFO];
  Int_t lptgen_mother[MAX_LPT_GENINFO];
  Int_t lptgen_befrad[MAX_LPT_GENINFO];
  Int_t lptgen_motherpdgid[MAX_LPT_GENINFO];
  Int_t lptgen_indinfo[MAX_LPT_GENINFO];
  Int_t lptgen_historycode[MAX_LPT_GENINFO];
  Float_t  lptgen_drmatch[MAX_LPT_GENINFO];
  Float_t lptgen_drmatchel[MAX_LPT_GENINFO];
  Float_t lptgen_drmatchmu[MAX_LPT_GENINFO];
  Float_t lptgen_drmatchph[MAX_LPT_GENINFO];
  Int_t lptgen_indrec[MAX_LPT_GENINFO];
  Int_t lptgen_indrecel[MAX_LPT_GENINFO];
  Int_t lptgen_indrecmu[MAX_LPT_GENINFO];
  Int_t lptgen_indrecph[MAX_LPT_GENINFO];

  TClonesArray *lptgeninfo_p4;
  TClonesArray *lptgen_p4;
  TClonesArray *lptgen_befrad_p4;

 private:
  //  double etCut_; 
  // edm::InputTag generatorColl;
  int debug_level;
  //  const HepMC::GenEvent* myGenEvent;
  double dR_min_for_matching; 

};

#endif
