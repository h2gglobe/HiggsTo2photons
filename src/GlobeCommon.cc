#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCommon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <iostream>

GlobeCommon::GlobeCommon(const edm::ParameterSet& iConfig) {
  generatorColl = iConfig.getParameter<edm::InputTag>("GeneratorColl");
}

void GlobeCommon::defineBranch(TTree* tree) {

  tree->Branch("event", &event, "event/I");
  tree->Branch("run", &run, "run/I");
  tree->Branch("process_id", &process_id, "process_id/I");
  tree->Branch("lumis", &lumis, "lumis/I");
  tree->Branch("bx", &bx, "bx/I");
  tree->Branch("pthat", &pthat, "pthat/F");
  tree->Branch("weight", &weight, "weight/F");
}

bool GlobeCommon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  static int nerror=0;

  edm::Handle<GenEventInfoProduct> HepMCEvt;
  edm::Handle<double> weightHandle;
  event = iEvent.id().event();
  run = iEvent.id().run();  
  lumis = iEvent.luminosityBlock();
  bx = iEvent.bunchCrossing();

  // add event PTHAT
  iEvent.getByLabel(generatorColl, HepMCEvt);
  iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
  
  if (HepMCEvt.isValid()) {

    const GenEventInfoProduct* event = HepMCEvt.product();
    pthat = event->qScale();
    process_id = event->signalProcessID();

    if (event->weights().size() > 1) {
      if (nerror++<10) {
	      std::cout << "GlobeCommon: more than 1 weight." << std::endl;
      }
    }

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
