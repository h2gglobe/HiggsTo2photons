#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeAnalyzer.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

GlobeAnalyzer::GlobeAnalyzer(const edm::ParameterSet& iConfig) {
  
  //jobmaker = new std::string("");

  fileName = iConfig.getParameter<std::string>("RootFileName");
  jobmaker = iConfig.getParameter<std::string>("JobMaker");
  
  doElectronStd = iConfig.getParameter<bool>("doElectron_std");

  doMuon = iConfig.getParameter<bool>("doMuon");

  doJetAlgo1 = iConfig.getParameter<bool>("doJet_algo1");
  doJetAlgo2 = iConfig.getParameter<bool>("doJet_algo2");
  doJetAlgo3 = iConfig.getParameter<bool>("doJet_algo3");
  doJetAlgoPF1 = iConfig.getParameter<bool>("doJet_algoPF1");
  doJetAlgoPF2 = iConfig.getParameter<bool>("doJet_algoPF2");
  doJetAlgoPF3 = iConfig.getParameter<bool>("doJet_algoPF3");

  doGenJetAlgo1 = iConfig.getParameter<bool>("doGenJet_algo1");
  doGenJetAlgo2 = iConfig.getParameter<bool>("doGenJet_algo2");
  doGenJetAlgo3 = iConfig.getParameter<bool>("doGenJet_algo3");

  doGenerator = iConfig.getParameter<bool>("doGenerator");
  doGenParticles = iConfig.getParameter<bool>("doGenParticles");
  doGenVertices = iConfig.getParameter<bool>("doGenVertices");

  doPhoton = iConfig.getParameter<bool>("doPhoton"); 
  doAllConversions = iConfig.getParameter<bool>("doAllConversions"); 
  doCaloTower = iConfig.getParameter<bool>("doCaloTower"); 
  doHcal = iConfig.getParameter<bool>("doHcal"); 
  doEcal = iConfig.getParameter<bool>("doEcal"); 

  doL1 = iConfig.getParameter<bool>("doL1");
  doHLT = iConfig.getParameter<bool>("doHLT");

  doVertices_std = iConfig.getParameter<bool>("doVertices_std"); 
  doVertices_nobs = iConfig.getParameter<bool>("doVertices_nobs"); 
  doMet = iConfig.getParameter<bool>("doMet"); 
  dotcMet = iConfig.getParameter<bool>("dotcMet");
  doPFMet = iConfig.getParameter<bool>("doPFMet");
  doSimHits = iConfig.getParameter<bool>("doSimHits"); 
  doSimTracks = iConfig.getParameter<bool>("doSimTracks"); 
  doTracks = iConfig.getParameter<bool>("doTracks"); 
  doGsfTracks = iConfig.getParameter<bool>("doGsfTracks"); 
  doTrackingParticles = iConfig.getParameter<bool>("doTrackingParticles"); 

  doEcalRecHits = iConfig.getParameter<bool>("doEcalRecHits");
  doReducedGen = iConfig.getParameter<bool>("doReducedGen");

  doLeptons = iConfig.getParameter<bool>("doLeptons");
  doHt = iConfig.getParameter<bool>("doHt");

  doPFCandidates = iConfig.getParameter<bool>("doPFCandidates");
  //doPAT = iConfig.getParameter<bool>("doPAT");

  doRho = iConfig.getParameter<bool>("doRho");
  doPileup = iConfig.getParameter<bool>("doPileup");

  debug_level = iConfig.getParameter<int>("Debug_Level");
  

  common = new GlobeCommon(iConfig);

  if(doGenerator && doGenParticles) {
    std::cout << "doGenerator and doGenParticles cannot be true at the same time." << std::endl;
    throw;
  }
    
  //GENERATOR
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeGenerator"<<std::endl;
  if (doGenerator)
    gen = new GlobeGenerator(iConfig);

  if (doGenParticles)
    genP = new GlobeGenParticles(iConfig);
  
  if (doGenVertices)
    genV = new GlobeGenVertices(iConfig);


  //SIMHITS
  if (doSimHits)
    simhits = new GlobeSimHits(iConfig);

  //SIMTRACKS
  if (doSimTracks)
    simtracks = new GlobeSimTracks(iConfig);

  //LEVEL 1 TRIGGER
  if (doL1)
    level1   = new GlobeL1(iConfig); 

  //HLTRIGGER
  if (doHLT)
    hlt = new GlobeHLT(iConfig); 
   
  //ECAL HCAL TRACKS VERTICES
  if (doEcal)
    ecalclusters = new GlobeEcalClusters(iConfig);

  if (doEcalRecHits)
    ecalrechits = new GlobeEcalHits(iConfig);

  if (doHcal)
    hcalhits = new GlobeHcal(iConfig);

  if (doCaloTower)
    calotowers = new GlobeCaloTowers(iConfig);

  if (doTracks)
    tracks = new GlobeTracks(iConfig);

  if (doGsfTracks)
    gsfTracks = new GlobeGsfTracks(iConfig);

  if (doTrackingParticles)
    trackingParticles = new GlobeTrackingParticles(iConfig);

  if (doVertices_std) 
    vertex_std   = new GlobeVertex(iConfig, "std");

  if (doVertices_nobs) 
    vertex_nobs   = new GlobeVertex(iConfig, "nobs");

  //PHOTONS
  if (doPhoton)
    photons = new GlobePhotons(iConfig);

  // Conversions
  if (doAllConversions)
    allConversions = new GlobeConversions(iConfig);

  //ELECTRONS

  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeElectrons"<<std::endl;
  if (doElectronStd)
    std_electrons = new GlobeElectrons(iConfig, "std");
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeElectrons"<<std::endl;

  //MUONS
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeMuons"<<std::endl;
  if (doMuon)
    muons = new GlobeMuons(iConfig);

  //MET
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeMET"<<std::endl;
  if (doMet)
    met = new GlobeMET(iConfig, "met"); 
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeTCMET"<<std::endl;
  if (dotcMet)
    tcmet = new GlobeMET(iConfig, "tcmet");
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobePFMET"<<std::endl;
  if (doPFMet)
    pfmet = new GlobeMET(iConfig, "pfmet");

  //GENJETS
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: coutall GlobeGenJets"<<std::endl;
  if (doGenJetAlgo1)
    algo1_genJets = new GlobeGenJets(iConfig, "algo1");
  if (doGenJetAlgo2)
    algo2_genJets = new GlobeGenJets(iConfig, "algo2");
  if (doGenJetAlgo3)
    algo3_genJets = new GlobeGenJets(iConfig, "algo3");

  //JETS
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeJets"<<std::endl;
  if (doJetAlgo1)
    algo1_jets = new GlobeJets(iConfig, "algo1");
  if (doJetAlgo2)
    algo2_jets = new GlobeJets(iConfig, "algo2");
  if (doJetAlgo3)
    algo3_jets = new GlobeJets(iConfig, "algo3");
  if (doJetAlgoPF1)
    algoPF1_jets = new GlobeJets(iConfig, "algoPF1");
  if (doJetAlgoPF2)
    algoPF2_jets = new GlobeJets(iConfig, "algoPF2");
  if (doJetAlgoPF3)
    algoPF3_jets = new GlobeJets(iConfig, "algoPF3");
  
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeSelector"<<std::endl;
  selector = new GlobeSelector(iConfig);

  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeLeptons"<<std::endl;
  if (doLeptons)
    leptons = new GlobeLeptons();
  
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeHt"<<std::endl;
  if (doHt)
    ht = new GlobeHT(iConfig);

  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeReducedGen"<<std::endl;
  if (doReducedGen)
    reducedgen = new GlobeReducedGen(iConfig);

  if (doPFCandidates)
    pfCandidates = new GlobePFCandidates(iConfig);

  if (doRho)
    rho = new GlobeRho(iConfig);

  if (doPileup)
    pileup = new GlobePileup(iConfig);

  //if (doPAT)
  //  pat = new GlobePAT(iConfig);

  if(debug_level > 9) std::cout<<"GlobeAnalyzer: readConfiguration"<<std::endl;
  readConfiguration(iConfig);
}

GlobeAnalyzer::~GlobeAnalyzer() {}

void GlobeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

 if(debug_level > 9) std::cout << "GlobeAnalyzer: Begin" << std::endl;
  
  tot_events++;

  common->analyze(iEvent, iSetup);

  //PHOTONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: photons" << std::endl;
  if (doPhoton)
    photons->analyze(iEvent, iSetup);

  //ALL CONVERSIONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: conversions" << std::endl;
  if (doAllConversions)
    allConversions->analyze(iEvent, iSetup);



  //ELECTRONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: std_electrons" << std::endl;
  if (doElectronStd)
    std_electrons->analyze(iEvent, iSetup);

  //MUONS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: muons" << std::endl;
  if (doMuon)
    muons->analyze(iEvent, iSetup);
  
  //Leptons
  if (doLeptons)
    leptons->Zero();

  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: leptons->addmuons" << std::endl;
  if(doMuon && doLeptons)
    leptons->addMuons(muons);

  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: leptons->addelectrons" << std::endl;
  if(doElectronStd && doLeptons)
    leptons->addElectrons(std_electrons);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: leptons->addphotons" << std::endl;
  if(doPhoton && doLeptons)
    leptons->addPhotons(photons);
 


  //GENERATOR
  if(debug_level > 2) std::cout << "GlobeAnalyzer: gen" << std::endl;
  if (doGenerator){
    gen->analyze(iEvent, iSetup);
  }

  //GENPARTICLES
  if(debug_level > 2) std::cout << "GlobeAnalyzer: genP" << std::endl;
  if (doGenParticles){
    genP->analyze(iEvent, iSetup);
  }

  //GENVERTICES 
  if(debug_level > 2) std::cout << "GlobeAnalyzer: genV" << std::endl;
  if (doGenVertices){
    genV->analyze(iEvent, iSetup);
  }

  //SIMHITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: simhits" << std::endl;
  if (doSimHits)
    simhits->analyze(iEvent, iSetup);

  //SIMTRACKS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: simtracks" << std::endl;
  if (doSimTracks)
    simtracks->analyze(iEvent, iSetup);

  //LEVEL 1 TRIGGER
  if(debug_level > 2) std::cout << "GlobeAnalyzer: level1" << std::endl;
  if (doL1)
    level1->analyze(iEvent, iSetup);

  //HLTRIGGER
  if(debug_level > 2) std::cout << "GlobeAnalyzer: hlt" << std::endl;
  if (doHLT)
    hlt->analyze(iEvent, iSetup);

  //TRACKS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: tracks" << std::endl;
  if (doTracks)
    tracks->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: gsftracks" << std::endl;
  if (doGsfTracks)
    gsfTracks->analyze(iEvent, iSetup);

  //TRACKING PARTICLES
  if(debug_level > 2) std::cout << "GlobeAnalyzer: trackingparticles" << std::endl;
  if (doTrackingParticles) {
    trackingParticles->analyze(iEvent, iSetup,tracks);
    tracks->GetAssociatedTrackingParticleIndex(iEvent,iSetup,trackingParticles);
  }

  //ECAL REC HITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: ecalrechits" << std::endl;
  if (doEcalRecHits) {
    if (doElectronStd && doMuon && doPhoton && doLeptons)
      ecalrechits->analyze(iEvent, iSetup, leptons, std_electrons, muons, photons);
    else
      std::cout << "EcalRecHits needs Electrons, Muons, Leptons and Photons." << std::endl;
  }

  //ECAL CLUSTERS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: ecalclusters" << std::endl;
  if (doEcal)
    ecalclusters->analyze(iEvent, iSetup);

  //HCAL HITS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: hcalhits" << std::endl;
  if (doHcal) {
    if (doElectronStd && doMuon && doPhoton && doTracks && doLeptons) 
      hcalhits->analyze(iEvent, iSetup, leptons, std_electrons, muons, photons, tracks);
    else
      std::cout << "HcalRecHits needs Electrons, Muons, Photons, Leptons and Tracks." << std::endl;
  }

  //CALO TOWERS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: calotowers" << std::endl;
  if (doCaloTower)
    calotowers->analyze(iEvent, iSetup);

  //VERTEX
  if(debug_level > 2) std::cout << "GlobeAnalyzer: vertex_std" << std::endl;
  if (doVertices_std) 
    vertex_std->analyze(iEvent, iSetup);
    
  if(debug_level > 2) std::cout << "GlobeAnalyzer: vertex_nobs" << std::endl;
  if (doVertices_nobs) 
    vertex_nobs->analyze(iEvent, iSetup);


  //MET
  if(debug_level > 2) std::cout << "GlobeAnalyzer: met" << std::endl;
  if (doMet)
    met->analyze(iEvent, iSetup);
  if (dotcMet)
    tcmet->analyze(iEvent, iSetup);
  if (doPFMet)
    pfmet->analyze(iEvent, iSetup);

  //GENJETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo1genjet" << std::endl;
  if (doGenJetAlgo1)
    algo1_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo2genjet" << std::endl;
  if (doGenJetAlgo2)
    algo2_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo3genjet" << std::endl;
  if (doGenJetAlgo3)
    algo3_genJets->analyze(iEvent, iSetup);

  //JETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo1_jets" << std::endl;
  if (doJetAlgo1)
    algo1_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo2_jets" << std::endl;
  if (doJetAlgo2)
    algo2_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algo3_jets" << std::endl;
  if (doJetAlgo3)
    algo3_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf1_jets" << std::endl;
  if (doJetAlgoPF1)
    algoPF1_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf2_jets" << std::endl;
  if (doJetAlgoPF2)
    algoPF2_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: algopf3_jets" << std::endl;
  if (doJetAlgoPF3)
    algoPF3_jets->analyze(iEvent, iSetup);

  //HT
  if (doHt) {
    bool doLeptonHT=doJetAlgo1 && doMet && doPhoton && doLeptons;
    if(debug_level > 2) std::cout << "GlobeAnalyzer: leptonHT" << std::endl;
    if(doLeptonHT)
      ht->fillLeptonHT(algo1_jets, met, leptons);
    
    bool doCaloTowerHT=doJetAlgo1 && doMet && doPhoton && doLeptons;
    if(debug_level > 2) std::cout << "GlobeAnalyzer: caloHT" << std::endl;
    if(doCaloTowerHT)
      ht->fillCaloTowerHT(met, calotowers);

    if(debug_level > 2) std::cout << "GlobeAnalyzer: trackHT" << std::endl;
    if(doTracks)
      ht->fillTrackHT(iEvent);
  }

  if(doReducedGen && doGenerator && doLeptons) {
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: reducedgen" << std::endl;
    reducedgen->fillRedGenList(gen, leptons);
  }

  if(doReducedGen && doGenParticles && doLeptons) {
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: reducedgenparticles" << std::endl;
    reducedgen->fillRedGenList(genP, leptons);
  }
  
  // PF CANDIDATES
  if (doPFCandidates)
    pfCandidates->analyze(iEvent, iSetup, tracks, muons);

  //PAT
  //if (doPAT) {
  //  if(debug_level > 2) std::cout << "GlobeAnalyzer: PAT" << std::endl;
  //  pat->analyze(iEvent, iSetup, std_electrons, photons, algo1_jets);
  //}

  if (doRho)
    rho->analyze(iEvent, iSetup);

  if (doPileup)
    pileup->analyze(iEvent, iSetup);

  

  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: selectorbits" << std::endl;

  if (doGenParticles || doGenerator)
    selector_bits = selector->select(std_electrons, muons, photons, gen, leptons, reducedgen).to_ulong();
  else
    selector_bits = selector->select(std_electrons, muons, photons).to_ulong();
  
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: selectorbits = " << selector_bits << std::endl;

  if (selector_bits > 0) {
    sel_events++;
    
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: fill my tree!" << std::endl;

    // fill the tree
    tree->Fill();
  }

  if(debug_level > 9) std::cout << "GlobeAnalyzer: End" << std::endl;
}

void GlobeAnalyzer::beginJob() { 

  // define root file and root tree
  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event", "Event data");
  tree2 = new TTree("global_variables", "Global Parameters"); // need a different tree to fill once per job

  common->defineBranch(tree);
  
  if (doPhoton)
    photons->defineBranch(tree);

  if (doAllConversions)
    allConversions->defineBranch(tree);  

  if (doEcalRecHits)
    ecalrechits->defineBranch(tree);
  if (doEcal)
    ecalclusters->defineBranch(tree);
  if (doCaloTower)
    calotowers->defineBranch(tree);
  if (doHcal)
    hcalhits->defineBranch(tree);

  if (doL1)
    level1->defineBranch(tree);

  if (doHLT)
    hlt->defineBranch(tree);

  if (doVertices_std)
    vertex_std->defineBranch(tree);
  if (doVertices_nobs)
    vertex_nobs->defineBranch(tree);

  if (doMet)
    met->defineBranch(tree);
  if(dotcMet)
    tcmet->defineBranch(tree);
  if(doPFMet)
    pfmet->defineBranch(tree);

  if (doSimHits)
    simhits->defineBranch(tree);
  if (doSimTracks)
    simtracks->defineBranch(tree);

  if (doTracks)
    tracks->defineBranch(tree);
  if (doGsfTracks)
    gsfTracks->defineBranch(tree);

  if (doTrackingParticles)
    trackingParticles->defineBranch(tree);
 
  // define object specific branches
  // ELECTRONS
  if (doElectronStd)
    std_electrons->defineBranch(tree);
 
  // MUONS
  if (doMuon)
    muons->defineBranch(tree);

  // JETS
  if (doJetAlgo1)
    algo1_jets->defineBranch(tree);
  if (doJetAlgo2)
    algo2_jets->defineBranch(tree);
  if (doJetAlgo3)
    algo3_jets->defineBranch(tree);
  if (doJetAlgoPF1)
    algoPF1_jets->defineBranch(tree);
  if (doJetAlgoPF2)
    algoPF2_jets->defineBranch(tree);
  if (doJetAlgoPF3)
    algoPF3_jets->defineBranch(tree);

  // GEN
  if (doGenerator)
    gen->defineBranch(tree);

  if (doGenParticles)
    genP->defineBranch(tree);

  if (doGenVertices)
    genV->defineBranch(tree);

  
  // GEN JETS
  if (doGenJetAlgo1)
    algo1_genJets->defineBranch(tree);
  if (doGenJetAlgo2)
    algo2_genJets->defineBranch(tree);
  if (doGenJetAlgo3)
    algo3_genJets->defineBranch(tree);

  if (doLeptons)
    leptons->defineBranch(tree);
  
  if (doHt)
    ht->defineBranch(tree);
  
  // REDUCED GEN LIST
  if (doReducedGen) 
    reducedgen->defineBranch(tree);

  if (doPFCandidates)
     pfCandidates->defineBranch(tree);
  
  //if (doPAT)
  //  pat->defineBranch(tree);


  if (doRho)
     rho->defineBranch(tree);

  if (doPileup)
     pileup->defineBranch(tree);
  
  defineBranch();
  
  tot_events = 0;
  sel_events = 0;
}

void GlobeAnalyzer::endJob() { 

  fillTree();
  file->Write();
  file->Close();
}

void GlobeAnalyzer::readConfiguration(const edm::ParameterSet& iConfig) {
  parameters = new std::vector<std::string>;
  parameters->push_back(iConfig.dump());
}

void GlobeAnalyzer::defineBranch() {

  tree->Branch("selector_bits", &selector_bits, "selector_bits/I");

  tree2->Branch("version", &version, "version/I");
  tree2->Branch("type", &type, "type/I");
  tree2->Branch("tot_events", &tot_events, "tot_events/I");
  tree2->Branch("sel_events", &sel_events, "sel_events/I");
  tree2->Branch("parameters", "std::vector<std::string>", &parameters); 
  tree2->Branch("jobmaker", "std::string", &jobmaker); 
}

void GlobeAnalyzer::fillTree() {

  version = H2G_VERSION;
  type = 0;

  //if (getenv("JOBMAKER") != NULL) 
  //  jobmaker = new std::string(getenv("JOBMAKER"));
  //else 
  //  jobmaker = new std::string("jobmaker unknown");

  tree2->Fill();
}

DEFINE_FWK_MODULE(GlobeAnalyzer);
