// Package:   multiSCfilter 
// Class:     multiSCfilter 
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Jun 22 13:56:58 CEST 2010
// copied/altered version of InvariantMassFilter.cc

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"

class multiSCfilter : public edm::EDFilter {
public:
  explicit multiSCfilter(const edm::ParameterSet&);
  ~multiSCfilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag scColl_;
  
  Int_t nSCrequire_;
  Float_t minEt_;
  Float_t maxEta_;
};


multiSCfilter::multiSCfilter(const edm::ParameterSet& iConfig) {
  
  scColl_  = iConfig.getParameter<edm::InputTag>("superClusterCollection");

  nSCrequire_ =  iConfig.getParameter<int>("nSCrequire");
  minEt_      =  iConfig.getParameter<double>("minEt");
  maxEta_     =  iConfig.getParameter<double>("maxEta");
}

multiSCfilter::~multiSCfilter()
{}

bool multiSCfilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //bool sel = false;
  
  edm::Handle<reco::SuperClusterCollection> scH;
  iEvent.getByLabel(scColl_, scH);

  int nSCpassing=0;

  for(reco::SuperClusterCollection::const_iterator sc = scH->begin(); sc != scH->end(); sc++) {
    math::PtEtaPhiMLorentzVector scp4(sin(sc->position().theta())*sc->energy(), sc->position().eta(), sc->position().phi(), 0);
    if (scp4.Et() > minEt_  && fabs(scp4.Eta()) < maxEta_)nSCpassing++;
    if(nSCpassing >= nSCrequire_) return true;
  }

  return false;
}

void multiSCfilter::beginJob()
{}

void multiSCfilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(multiSCfilter);
