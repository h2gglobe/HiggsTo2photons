#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeAnalyzer.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

GlobeAnalyzer::GlobeAnalyzer(const edm::ParameterSet& iConfig) {
  
  fileName = iConfig.getParameter<std::string>("RootFileName");
  
  doElectronStd = iConfig.getParameter<bool>("doElectron_std");

  doMuon = iConfig.getParameter<bool>("doMuon");

  doJetIt5 = iConfig.getParameter<bool>("doJet_it5");
  doJetIt7 = iConfig.getParameter<bool>("doJet_it7");
  doJetMid = iConfig.getParameter<bool>("doJet_mid");
  doJetIt5PF = iConfig.getParameter<bool>("doJet_it5pf");
  doJetSis5PF = iConfig.getParameter<bool>("doJet_sis5pf");
  doJetKt4PF = iConfig.getParameter<bool>("doJet_kt4pf");

  doGenJetIt5 = iConfig.getParameter<bool>("doGenJet_it5");
  doGenJetIt7 = iConfig.getParameter<bool>("doGenJet_it7");
  doGenJetMid = iConfig.getParameter<bool>("doGenJet_mid");

  doGenerator = iConfig.getParameter<bool>("doGenerator");
  doGenParticles = iConfig.getParameter<bool>("doGenParticles");

  doPhoton = iConfig.getParameter<bool>("doPhoton"); 
  doCaloTower = iConfig.getParameter<bool>("doCaloTower"); 
  doHcal = iConfig.getParameter<bool>("doHcal"); 
  doEcal = iConfig.getParameter<bool>("doEcal"); 

  doL1 = iConfig.getParameter<bool>("doL1");
  doHLT = iConfig.getParameter<bool>("doHLT");

  doVertices_std = iConfig.getParameter<bool>("doVertices_std"); 
  doVertices_pix = iConfig.getParameter<bool>("doVertices_pix"); 
  doVtxCompat = iConfig.getParameter<bool>("doVtxCompat"); 
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

  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  fullHLT = iConfig.getParameter<edm::ParameterSet>("HLTParameters").getParameter<bool>("FullHLT");
  if(fullHLT){
    theElHLTLabels  = iConfig.getParameter<std::vector<edm::InputTag> >("ElectronHLTLabels");
    theMuHLTLabels  = iConfig.getParameter<std::vector<edm::InputTag> >("MuonHLTLabels");
    thePhHLTLabels  = iConfig.getParameter<std::vector<edm::InputTag> >("PhotonHLTLabels");
  }
  
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
  if (doVertices_pix) 
    vertex_pix   = new GlobeVertex(iConfig, "pix");

  if (doVtxCompat) 
    vtxcompat   = new GlobeVtxCompat(iConfig);

  //PHOTONS
  if (doPhoton)
    photons = new GlobePhotons(iConfig);

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
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeGenJets"<<std::endl;
  if (doGenJetIt5)
    it5_genJets = new GlobeGenJets(iConfig, "it5");
  if (doGenJetIt7)
    it7_genJets = new GlobeGenJets(iConfig, "it7");
  if (doGenJetMid)
    mid_genJets = new GlobeGenJets(iConfig, "mid");

  //JETS
  if(debug_level > 9) std::cout<<"GlobeAnalyzer: call GlobeJets"<<std::endl;
  if (doJetIt5)
    it5_jets = new GlobeJets(iConfig, "it5");
  if (doJetIt7)
    it7_jets = new GlobeJets(iConfig, "it7");
  if (doJetMid)
    mid_jets = new GlobeJets(iConfig, "mid");
  if (doJetIt5PF)
    it5pf_jets = new GlobeJets(iConfig, "it5pf");
  if (doJetSis5PF)
    sis5pf_jets = new GlobeJets(iConfig, "sis5pf");
  if (doJetKt4PF)
    kt4pf_jets = new GlobeJets(iConfig, "kt4pf");
  
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
    
  if(debug_level > 2) std::cout << "GlobeAnalyzer: vertex_pix" << std::endl;
  if (doVertices_pix) 
    vertex_pix->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: vtxcompat" << std::endl;
  if (doVtxCompat) { 
    if (doElectronStd && doMuon && doPhoton && doTracks && doLeptons) 
      vtxcompat->analyze(iEvent, iSetup, leptons, std_electrons, muons, tracks);
    else
      std::cout << "VtxCompat needs Electrons, Muons, Photons, Leptons and Tracks." << std::endl;
  }  

  //MET
  if(debug_level > 2) std::cout << "GlobeAnalyzer: met" << std::endl;
  if (doMet)
    met->analyze(iEvent, iSetup);
  if (dotcMet)
    tcmet->analyze(iEvent, iSetup);
  if (doPFMet)
    pfmet->analyze(iEvent, iSetup);

  //GENJETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: it5genjet" << std::endl;
  if (doGenJetIt5)
    it5_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: it7genjet" << std::endl;
  if (doGenJetIt7)
    it7_genJets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: midgenJets" << std::endl;
  if (doGenJetMid)
    mid_genJets->analyze(iEvent, iSetup);

  //JETS
  if(debug_level > 2) std::cout << "GlobeAnalyzer: t5_jets" << std::endl;
  if (doJetIt5)
    it5_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: t7_jets" << std::endl;
  if (doJetIt7)
    it7_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: mid_jets" << std::endl;
  if (doJetMid)
    mid_jets->analyze(iEvent, iSetup);

  if(debug_level > 2) std::cout << "GlobeAnalyzer: pfit5_jets" << std::endl;
  if (doJetIt5PF)
    it5pf_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: pfsis5_jets" << std::endl;
  if (doJetSis5PF)
    sis5pf_jets->analyze(iEvent, iSetup);
  if(debug_level > 2) std::cout << "GlobeAnalyzer: pfkt4_jets" << std::endl;
  if (doJetKt4PF)
    kt4pf_jets->analyze(iEvent, iSetup);

  //HT
  if (doHt) {
    bool doLeptonHT=doJetIt5 && doMet && doPhoton && doLeptons;
    if(debug_level > 2) std::cout << "GlobeAnalyzer: leptonHT" << std::endl;
    if(doLeptonHT)
      ht->fillLeptonHT(it5_jets, met, leptons);
    
    bool doCaloTowerHT=doJetIt5 && doMet && doPhoton && doLeptons;
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
  
  if(debug_level > 2) 
    std::cout << "GlobeAnalyzer: selectorbits" << std::endl;

  if (doGenParticles || doGenerator)
    selector_bits = selector->select(std_electrons, muons, photons, gen, leptons, reducedgen).to_ulong();
  else
    selector_bits = selector->select(std_electrons, muons, photons).to_ulong();
  
  if (selector_bits > 0) {
    sel_events++;
    
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
  if (doVertices_pix)
    vertex_pix->defineBranch(tree);

  if (doVtxCompat)
    vtxcompat->defineBranch(tree);

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
  if (doJetIt5)
    it5_jets->defineBranch(tree);
  if (doJetIt7)
    it7_jets->defineBranch(tree);
  if (doJetMid)
    mid_jets->defineBranch(tree);
  if (doJetIt5PF)
    it5pf_jets->defineBranch(tree);
  if (doJetSis5PF)
    sis5pf_jets->defineBranch(tree);
  if (doJetKt4PF)
    kt4pf_jets->defineBranch(tree);

  // GEN
  if (doGenerator)
    gen->defineBranch(tree);

  if (doGenParticles)
    genP->defineBranch(tree);

  
  // GEN JETS
  if (doGenJetIt5)
    it5_genJets->defineBranch(tree);
  if (doGenJetIt7)
    it7_genJets->defineBranch(tree);
  if (doGenJetMid)
    mid_genJets->defineBranch(tree);

  if (doLeptons)
    leptons->defineBranch(tree);
  
  if (doHt)
    ht->defineBranch(tree);
  
  // REDUCED GEN LIST
  if (doReducedGen) 
    reducedgen->defineBranch(tree);
  
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
  if(fullHLT){
    hlt_path_names = new std::vector<std::string>;
    tree2->Branch("hlt_path_names", "std::vector<std::string>", &hlt_path_names);
  }
}

void GlobeAnalyzer::fillTree() {

  version = H2G_VERSION;
  type = 0;
  
  if(fullHLT){
    for(unsigned int i=0; i<theElHLTLabels.size(); i++) 
      hlt_path_names->push_back(theElHLTLabels[i].label());
    for(unsigned int i=0; i<theMuHLTLabels.size(); i++) 
      hlt_path_names->push_back(theMuHLTLabels[i].label());
    for(unsigned int i=0; i<thePhHLTLabels.size(); i++) 
      hlt_path_names->push_back(thePhHLTLabels[i].label());
  }
  
  tree2->Fill();
}

DEFINE_FWK_MODULE(GlobeAnalyzer);
