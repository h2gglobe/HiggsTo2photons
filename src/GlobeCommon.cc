#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCommon.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <iostream>

GlobeCommon::GlobeCommon(const edm::ParameterSet& iConfig) {
  generatorColl = iConfig.getParameter<edm::InputTag>("GeneratorColl");
  doParticleGun =iConfig.getParameter<bool>("doParticleGun");
}

void GlobeCommon::defineBranch(GlobeAnalyzer* ana) {

  ana->Branch("event", &event, "event/I");
  ana->Branch("run", &run, "run/I");
  ana->Branch("process_id", &process_id, "process_id/I");
  ana->Branch("lumis", &lumis, "lumis/I");
  ana->Branch("bx", &bx, "bx/I");
  ana->Branch("pthat", &pthat, "pthat/F");
  ana->Branch("weight", &weight, "weight/F");
}

void GlobeCommon::defineLumiBranch(TTree* ana) {

  ana->Branch("run", &run, "run/I");
  ana->Branch("lumis", &lumis, "lumis/I");
}

bool GlobeCommon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  static int nerror=0;

  edm::Handle<GenEventInfoProduct> HepMCEvt;
  edm::Handle<double> weightHandle;
  event = iEvent.id().event();
  run = iEvent.id().run();  
  lumis = iEvent.luminosityBlock();
  bx = iEvent.bunchCrossing();
  //std::cout<<"event run "<<event<<" "<<run<<std::endl;
  // add event PTHAT
  iEvent.getByLabel(generatorColl, HepMCEvt);
  iEvent.getByLabel ("generator", "weight", weightHandle);
  
  if (HepMCEvt.isValid()) {

    const GenEventInfoProduct* event = HepMCEvt.product();
    pthat = event->qScale();
    
    if (doParticleGun)
      process_id = -1;
    else
      process_id = event->signalProcessID();

    if (event->weights().size() > 1) {
      if (nerror++<10) {
	      std::cout << "GlobeCommon: more than 1 weight." << std::endl;
      }
    }
    
    if (doParticleGun)
      weight = -1; 
    else
      weight = event->weights().front();
  } else {
    pthat = -1;
    process_id = -1;
    weight = 1;
  }
  
  if (weightHandle.isValid())  
    weight = (*weightHandle);

  return true;
}

void GlobeCommon::endLumiBlock(const edm::LuminosityBlock & l, const edm::EventSetup & es) {
  lumis = l.id().luminosityBlock();
  run   = l.id().run();
}
