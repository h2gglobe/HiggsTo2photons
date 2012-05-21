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
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class InvariantMassFilter : public edm::EDFilter {
public:
  explicit InvariantMassFilter(const edm::ParameterSet&);
  ~InvariantMassFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag elColl_;
  edm::InputTag scColl_;
  bool invMass_ ;
  Float_t minMassCut_;
  Float_t maxMassCut_;
};


InvariantMassFilter::InvariantMassFilter(const edm::ParameterSet& iConfig) {
  
  elColl_  = iConfig.getParameter<edm::InputTag>("electronCollection");
  scColl_  = iConfig.getParameter<edm::InputTag>("superClusterCollection");
  invMass_ = iConfig.getParameter<bool> ("ElectronElectronInvMass");
  
  minMassCut_ =  iConfig.getParameter<double>("minMassCut");
  maxMassCut_ =  iConfig.getParameter<double>("maxMassCut");
}

InvariantMassFilter::~InvariantMassFilter()
{}

bool InvariantMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<reco::SuperClusterCollection> scH;
  if (!invMass_)
    iEvent.getByLabel(scColl_, scH);

  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(elColl_, elH);

  if (!invMass_) {
    for(reco::GsfElectronCollection::const_iterator igsf = elH->begin(); igsf != elH->end(); igsf++) {
      math::XYZTLorentzVector l1 = igsf->p4();
      for(reco::SuperClusterCollection::const_iterator sc = scH->begin(); sc != scH->end(); sc++) {
	math::PtEtaPhiMLorentzVector l2(sin(sc->position().theta())*sc->energy(),
					sc->position().eta(),
					sc->position().phi(),
					0);
	if (((l1+l2).M() < maxMassCut_) && ((l1+l2).M() > minMassCut_))
	  return true;
      }
    }
  } else {
    for(reco::GsfElectronCollection::const_iterator igsf1 = elH->begin(); igsf1 != elH->end(); igsf1++) {
      math::XYZTLorentzVector l1 = igsf1->p4();
      for(reco::GsfElectronCollection::const_iterator igsf2 = elH->begin(); igsf2 != elH->end(); igsf2++) {
	if (igsf1 == igsf2)
	  continue;
	math::XYZTLorentzVector l2 = igsf2->p4();
	if (((l1+l2).M() < maxMassCut_) && ((l1+l2).M() > minMassCut_))
	  return true;
      }
    }
  }
  
  return false;
}

void InvariantMassFilter::beginJob()
{}

void InvariantMassFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(InvariantMassFilter);
