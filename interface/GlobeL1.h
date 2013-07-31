#ifndef GLOBEL1_H
#define GLOBEL1_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include <iostream>

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeAnalyzer;

class GlobeL1 {
 public:
  
  GlobeL1(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeL1() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

  // variables
  Int_t l1emiso_n;
  //int l1emiso_rank[MAX_L1];
  Float_t l1emiso_et[MAX_L1]; 
  Float_t l1emiso_eta[MAX_L1]; 
  Float_t l1emiso_phi[MAX_L1];
  
  Int_t l1emnoniso_n;
  //int l1emnoniso_rank[MAX_L1];
  Float_t l1emnoniso_et[MAX_L1];
  Float_t l1emnoniso_eta[MAX_L1];
  Float_t l1emnoniso_phi[MAX_L1];

  Int_t l1cenjet_n;
  //int l1cenjet_rank[MAX_L1];
  Float_t l1cenjet_et[MAX_L1];
  Float_t l1cenjet_eta[MAX_L1];
  Float_t l1cenjet_phi[MAX_L1];
  
  Int_t l1forjet_n;
  //int l1forjet_rank[MAX_L1];
  Float_t l1forjet_et[MAX_L1];
  Float_t l1forjet_eta[MAX_L1];
  Float_t l1forjet_phi[MAX_L1];
  
  Int_t l1taujet_n;
  //int l1taujet_rank[MAX_L1];
  Float_t l1taujet_et[MAX_L1];
  Float_t l1taujet_eta[MAX_L1];
  Float_t l1taujet_phi[MAX_L1];
  
  //unsigned int l1met_n;
  //int l1met_rank[MAX_L1];
  Float_t l1met_et;
  //float l1met_etTotal; //CHECK UNUSED VARIABLE do you want to put it?
  Float_t l1met_phi;
  //float l1met_etHad; //CHECK UNUSED VARIABLE do you want to put it?

  Int_t l1mu_n;
  //int l1mu_rank[MAX_L1];
  Float_t l1mu_et[MAX_L1];
  Float_t l1mu_eta[MAX_L1];
  Float_t l1mu_phi[MAX_L1];
  
  std::vector<int>* l1bits_phy;
  std::vector<int>* l1bits_tec;
  std::map<std::string, int>* l1_labels;

 private:
  const char* nome;
  edm::InputTag l1EMIso;
  edm::InputTag l1EMNonIso;
  edm::InputTag l1CenJet;
  edm::InputTag l1ForJet;
  edm::InputTag l1TauJet;
  edm::InputTag l1EtMiss;
  edm::InputTag l1Mu;
  edm::InputTag m_gtReadoutRecord;
  edm::InputTag m_gtObjectMapRecord;

  int debug_level;
};

#endif
