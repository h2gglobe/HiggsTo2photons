// Package:    eIDFilter
// Class:      eIDFilter
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Wed Jun 16 19:04:58 CEST 2010

#include <memory>
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class eIDFilter : public edm::EDFilter {
public:
  explicit eIDFilter(const edm::ParameterSet&);
  ~eIDFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  bool useOr_;
  int selection_;
  std::vector<edm::InputTag> eIDCiCLabels_;
  std::vector<edm::InputTag> eIDVBTFLabels_;
  edm::InputTag elColl_;
};


eIDFilter::eIDFilter(const edm::ParameterSet& iConfig) {
  
  selection_ = iConfig.getParameter<int>("nElectrons");
  useOr_ = iConfig.getParameter<bool>("useOR");
  eIDCiCLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("eIDCiCLabels");
  eIDVBTFLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("eIDVBTFLabels");
  elColl_ = iConfig.getParameter<edm::InputTag>("electronCollection");

  produces<reco::GsfElectronCollection>();
}

eIDFilter::~eIDFilter()
{}

bool eIDFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::vector<edm::Handle<edm::ValueMap<float> > > eIDCiCs(eIDCiCLabels_.size()); 
  std::vector<edm::Handle<edm::ValueMap<float> > > eIDVBTFs(eIDVBTFLabels_.size()); 

  for (unsigned int i=0; i<eIDCiCLabels_.size(); i++)
    iEvent.getByLabel(eIDCiCLabels_[i], eIDCiCs[i]);

  for (unsigned int i=0; i<eIDVBTFLabels_.size(); i++)
    iEvent.getByLabel(eIDVBTFLabels_[i], eIDVBTFs[i]);

  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(elColl_, elH);


  std::auto_ptr<reco::GsfElectronCollection> electronCollection(new reco::GsfElectronCollection);

  int nSelectedEl = 0;
  for (unsigned int i=0; i<elH->size(); i++) {
    reco::GsfElectronRef electronRef(elH, i);

    bool selEl = false;
    if (!useOr_)
      selEl = true;

    for (unsigned int j=0; j<eIDCiCs.size(); j++) {
      int eID = (int)(*eIDCiCs[j])[electronRef];     
      if (useOr_)     
        selEl = selEl || (eID == 15);
      else
        selEl = selEl && (eID == 15);
    }

    for (unsigned int j=0; j<eIDVBTFs.size(); j++) {
      int eID = (int)(*eIDVBTFs[j])[electronRef];
      if (useOr_)
        selEl = selEl || (eID == 7);
      else
        selEl = selEl && (eID == 7);
    }

    if (selEl) {
      nSelectedEl++;
      electronCollection->push_back(*electronRef);
    }
  }

  iEvent.put(electronCollection);
  if (nSelectedEl >= selection_) {
    return true;
  }

  return false;
}


void eIDFilter::beginJob()
{}

void eIDFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(eIDFilter);
