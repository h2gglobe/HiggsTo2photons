// Package:    InvariantMassFilter
// Class:      InvariantMassFilter
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Jun 22 13:56:58 CEST 2010

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"

class PhotonInvariantMassFilter : public edm::EDFilter {
public:
  explicit PhotonInvariantMassFilter(const edm::ParameterSet&);
  ~PhotonInvariantMassFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag phoColl_;
  
  Float_t minMassCut_;
  Float_t maxMassCut_;
};


PhotonInvariantMassFilter::PhotonInvariantMassFilter(const edm::ParameterSet& iConfig) {
  
  phoColl_  = iConfig.getParameter<edm::InputTag>("photonCollection");
  
  minMassCut_ =  iConfig.getParameter<double>("minMassCut");
  maxMassCut_ =  iConfig.getParameter<double>("maxMassCut");
}

PhotonInvariantMassFilter::~PhotonInvariantMassFilter()
{}

bool PhotonInvariantMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(phoColl_, phoH);

  for(reco::PhotonCollection::const_iterator ipho = phoH->begin(); ipho+1 != phoH->end(); ++ipho ) {
	  for(reco::PhotonCollection::const_iterator jpho = ipho+1; jpho != phoH->end(); ++jpho ) {
		  float mass = (ipho->p4() + jpho->p4()).M();
		  if( mass > minMassCut_ && mass < maxMassCut_ ) { return true; }
	  }
  }
  
  return false;
}

void PhotonInvariantMassFilter::beginJob()
{}

void PhotonInvariantMassFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonInvariantMassFilter);
