#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBELEPTONS_H
#define GLOBELEPTONS_H

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <iostream>

class GlobeAnalyzer;

class GlobeLeptons {
  public:

    GlobeLeptons();
    virtual ~GlobeLeptons() {};

    void defineBranch(GlobeAnalyzer* ana);
    void fillList(GlobeElectrons * theElectrons, GlobeMuons * theMuons, GlobePhotons * thePhotons);
    void Zero();
    void addMuons(GlobeMuons * theMuons);
    void addElectrons(GlobeElectrons * theElectrons);
    void addPhotons(GlobePhotons * thePhotons);

    Int_t lpt_n;
    Int_t lpt_mu_n;
    Int_t lpt_el_n;
    Int_t lpt_emu_n;
    Int_t lpt_pho_n;
    Int_t lpt_pdgid[MAX_LEPTONS];
    Int_t lpt_ind[MAX_LEPTONS];
    Int_t lpt_duplicate[MAX_LEPTONS];
    Int_t lpt_indgen[MAX_LEPTONS];
    Float_t lpt_drmatch[MAX_LEPTONS];
    
    TClonesArray *lpt_p4;

  float leptons_cut_eta_rec;
  float leptons_cut_et_rec;

  private:
    int debug_level;
};

#endif
