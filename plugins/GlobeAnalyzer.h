// -*- C++ -*-
//
// Package:    GlobeAnalyzer
// Class:      GlobeAnalyzer
// 
/**\class GlobeAnalyzer GlobeAnalyzer.cc HiggsAnalysis/HiggsTo2photons/src/GlobeAnalyzer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Matteosan SANI
//         Created:  Thu Feb  7 10:14:43 CET 2008
// $Id: GlobeAnalyzer.h,v 1.2 2012/04/23 21:04:36 sani Exp $
//
//

#ifndef GLOBEANALYZER_H
#define GLOBEANALYZER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCommon.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeConversions.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCaloTowers.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHcal.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeL1.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeVertex.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMET.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePFCandidates.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimHits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGsfTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTrackingParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenerator.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenVertices.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenJets.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalHits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSelector.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHLT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeReducedGen.h"
//#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePAT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePdfWeights.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeRho.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePileup.h"

#include "TFile.h"
//#include "TRFIOFile.h"
#include "TTree.h"

// system include files
#include <memory>
#include <vector>
#include <string>

class GlobeAnalyzer : public edm::EDAnalyzer {
public:
  explicit GlobeAnalyzer(const edm::ParameterSet&);
  ~GlobeAnalyzer();
  void readConfiguration(const edm::ParameterSet& iConfig);
  void defineBranch();
  void fillTree();
  GlobeCommon* common;
  GlobePhotons* photons;
  GlobeConversions* allConversions;
  GlobePFCandidates* pfCandidates;
  GlobeEcalClusters* ecalclusters;
  GlobeMET* met, *tcmet, *pfmet; 
  GlobeCaloTowers* calotowers;
  GlobeHcal* hcalhits;
  GlobeL1* level1;
  GlobeVertex* vertex_std;
  GlobeVertex* vertex_nobs;
  GlobeSimHits* simhits;
  GlobeSimTracks* simtracks;
  GlobeTracks* tracks;
  GlobeGsfTracks* gsfTracks;
  GlobeTrackingParticles* trackingParticles;
  GlobeElectrons* std_electrons; //, *gge_electrons;
  GlobeMuons* global_muons, *tk_muons, *sta_muons, *muons;
  GlobeJets* algo1_jets, *algo2_jets, *algo3_jets, *algoPF1_jets, *algoPF2_jets, *algoPF3_jets ;
  GlobeGenerator* gen;
  GlobeGenParticles* genP;
  GlobeGenVertices* genV;
  GlobeGenJets* algo1_genJets, *algo2_genJets, *algo3_genJets;
  GlobeEcalHits* ecalrechits;
  GlobeHLT* hlt;
  GlobeSelector* selector;
  GlobeLeptons* leptons;
  GlobeHT* ht;
  //GlobePAT* pat;
  GlobeReducedGen* reducedgen;
  GlobeRho* rho1, *rho2, *rho3;
  GlobePileup* pileup;
  GlobePdfWeights* pdfweights;

private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  void endJob();

  std::string fileName;
      
  TFile *file;
  TTree *tree, *tree2, *lumitree;

  //std::vector<std::string> a, b;
  std::vector<std::string> *parameters;
  std::vector<std::string> *hlt_path_names, *reduced_path;
  std::vector<int>* reduced_index;
  std::string jobmaker;
  std::vector<std::string> globalCountersNames;

  std::vector<int> globalCounters, globalCountersPerLumi;
  
  int version, type, sel_events, tot_events; 
  int selector_bits;

  int debug_level;

  bool doElectronStd;
  bool doPFCandidates;
  bool doMuon;
  bool doMuonGlobal;
  bool doMuonTk;
  bool doMuonSta;
  bool doJetAlgo1;
  bool doJetAlgo2;
  bool doJetAlgo3;
  bool doJetAlgoPF1;
  bool doJetAlgoPF2;
  bool doJetAlgoPF3;
  bool doGenJetAlgo1;
  bool doGenJetAlgo2;
  bool doGenJetAlgo3;
  bool doGenerator;
  bool doGenParticles;
  bool doGenVertices;
  bool doCaloTower;
  bool doHcal;
  bool doEcal;
  bool doPhoton;
  bool doAllConversions;
  bool doL1;
  bool doVertices_std;
  bool doVertices_nobs;
  bool doMet;
  bool dotcMet;
  bool doPFMet;
  bool doSimHits;
  bool doSimTracks;
  bool doTracks;
  bool doGsfTracks;
  bool doTrackingParticles;
  bool doEcalRecHits;
  bool doTkRecHits;
  bool doHLT;
  bool doReducedGen;
  bool doLeptons;
  bool doHt;
  bool doPAT;
  bool doRho;
  bool doPileup;
  bool doPdfWeight;
  
  bool fullHLT;
  std::vector<edm::InputTag> theElHLTLabels;
  std::vector<edm::InputTag> theMuHLTLabels;
  std::vector<edm::InputTag> thePhHLTLabels;
  std::vector<edm::InputTag> theJetHLTLabels;

  Int_t nProcessedEvents;
};

#endif
