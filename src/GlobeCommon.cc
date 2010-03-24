#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCommon.h"

#ifdef CMSSW_VERSION_210
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#else
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#endif


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

  edm::Handle<edm::HepMCProduct> HepMCEvt;
  edm::Handle<double> weightHandle;
  event = iEvent.id().event();
  run = iEvent.id().run();  
  lumis = iEvent.luminosityBlock();
  bx = iEvent.bunchCrossing();

//#ifndef CMSSW_VERSION_210
//
//  edm::ProductID id_GenInfo, id_WeightInfo;
//  if (id_GenInfo = edm::findProductIDByLabel<edm::HepMCProduct>(iEvent, generatorColl), id_GenInfo.isValid()) {
//    // add event PTHAT
//    iEvent.getByLabel(generatorColl, HepMCEvt);
//  
//    const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
//    pthat = MCEvt->event_scale();  
//    process_id = MCEvt->signal_process_id();
//    if (MCEvt->weights().size() > 1) 
//      std::cout << "GlobeCommon: more than 1 weight." << std::endl;
//    weight = MCEvt->weights().front();
//  } else {
//    pthat = -1;  
//    process_id = -1;
//    weight = 1;
//  }
//  
//  if (id_WeightInfo = edm::findProductIDByLabel<double>(iEvent, "csa07EventWeightProducer","weight"), id_WeightInfo.isValid()) {
//    //iEvent.getByLabel(EventWeightProducer,EventWeightLabel,weightHandle);
//    iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
//    
//    weight = (*weightHandle);
//  }
//#else
  // add event PTHAT
  iEvent.getByLabel(generatorColl, HepMCEvt);
  iEvent.getByLabel ("csa07EventWeightProducer","weight", weightHandle);
  
  if (HepMCEvt.isValid()) {

    const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
    pthat = MCEvt->event_scale();  
    process_id = MCEvt->signal_process_id();

    if (MCEvt->weights().size() > 1) {
      if (nerror++<10) {
	std::cout << "GlobeCommon: more than 1 weight." << std::endl;
      }
    }

    weight = MCEvt->weights().front();
  } else {
    pthat = -1;
    process_id = -1;
    weight = 1;
  }
  
  if (weightHandle.isValid())  
    weight = (*weightHandle);

//#endif

  return true;
}
