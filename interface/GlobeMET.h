#ifndef GLOBEMET_H
#define GLOBEMET_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/METReco/interface/CaloMET.h" 

#include <iostream>

class GlobeAnalyzer;

class GlobeMET {
 public:
  
  enum CorrectionType {NoCaloCorrection, CrossedEnergyCorrection, S9EnergyCorrection, ExpectedMipEnergyCorrection};
  
  GlobeMET(const edm::ParameterSet&, const char* );
  virtual ~GlobeMET() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  void correctMETmuons(const edm::Event&, float&, float&, CorrectionType);
  void correctedJetMET(const edm::Event&, float&, float&, const float);

  // variables
  Float_t met_met;
  Float_t met_phi;
  Float_t met_met_nocalo;
  Float_t met_phi_nocalo;
  Float_t met_met_crossed;
  Float_t met_phi_crossed;
  Float_t met_met_s9;
  Float_t met_phi_s9;
  Float_t met_met_mip;
  Float_t met_phi_mip;
  Float_t met_met_jet;
  Float_t met_phi_jet;
  Float_t met_tcmet;
  Float_t met_phi_tcmet;
  Float_t met_pfmet;
  Float_t met_phi_pfmet;
  Float_t met_sumet_pfmet;
  Float_t met_mEtSig_pfmet;
  Float_t met_significance_pfmet;
  
  Float_t met_pfmetType1;
  Float_t met_phi_pfmetType1;
  Float_t met_sumet_pfmetType1;
  Float_t met_mEtSig_pfmetType1;
  Float_t met_significance_pfmetType1;

 private:
  const char* nome;
  edm::InputTag caloMETColl, muonGlobalColl, jetColl, tcMETColl, pfMETColl, pfType1METColl;
  int debug_level;
};

#endif
