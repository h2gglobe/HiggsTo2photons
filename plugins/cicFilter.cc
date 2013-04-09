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

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"

class cicFilter : public edm::EDFilter {
public:
  explicit cicFilter(const edm::ParameterSet&);
  ~cicFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  float SumPt2(reco::VertexRef vtx);
  int ChooseVertex(edm::Handle<reco::VertexCollection> vtxHandle);
    
  int selection_;
  bool useOR_;
  edm::InputTag phoColl_; 
  edm::InputTag pfColl_;
  edm::InputTag eleColl_;
  edm::InputTag vtxColl_;
  edm::InputTag tkColl_;
  edm::InputTag rhoColl_;

  int cutLevel_;

  CiCPhotonID* photonID;
};

cicFilter::cicFilter(const edm::ParameterSet& iConfig) {
  
  selection_ = iConfig.getParameter<int>("nPhotons");
  phoColl_   = iConfig.getParameter<edm::InputTag>("PhotonCollection");
  pfColl_    = iConfig.getParameter<edm::InputTag>("PFCollection");
  eleColl_   = iConfig.getParameter<edm::InputTag>("ElectronCollection");
  vtxColl_   = iConfig.getParameter<edm::InputTag>("VertexCollection");
  tkColl_    = iConfig.getParameter<edm::InputTag>("TrackCollection");
  rhoColl_   = iConfig.getParameter<edm::InputTag>("RhoCollection");

  useOR_     = iConfig.getParameter<bool>("useOR");
  cutLevel_  = iConfig.getParameter<int>("CutLevel");

  photonID = new CiCPhotonID(iConfig);

  produces<reco::PhotonCollection>();
}

cicFilter::~cicFilter()
{}

bool cicFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(phoColl_, phoH);

  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel(pfColl_, pfHandle);
  
  edm::Handle<reco::GsfElectronCollection> eleHandle;
  iEvent.getByLabel(eleColl_, eleHandle);

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByLabel(vtxColl_, vtxHandle);
  
  edm::Handle<reco::TrackCollection> tkHandle;
  iEvent.getByLabel(tkColl_, tkHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoColl_, rhoHandle);
  double rho = *(rhoHandle.product());

  photonID->configure(vtxHandle, tkHandle, eleHandle, pfHandle, rho);
  
  
  std::auto_ptr<reco::PhotonCollection> phoCollection(new reco::PhotonCollection);

  int nSelectedPho = 0;
  
  int chosenVertex = ChooseVertex(vtxHandle);

  if (chosenVertex == -1) 
    return false;

  for (unsigned int i=0; i<phoH->size(); i++) {

    reco::PhotonRef phoRef(phoH, i);
    bool pho_id_4cat   = photonID->PhotonID(4, phoRef, chosenVertex, CiCPhotonID::CiCPhotonIDLevel(cutLevel_));
    bool pho_id_6catpf = photonID->PhotonIDPF(6, phoRef, chosenVertex, CiCPhotonID::CiCPhotonIDLevel(cutLevel_));

    if (useOR_ and (pho_id_4cat or pho_id_6catpf)) {
      nSelectedPho++;
      phoCollection->push_back(*phoRef);
    }

    if (!useOR_ and (pho_id_4cat and pho_id_6catpf)) {
      nSelectedPho++;
      phoCollection->push_back(*phoRef);
    }
  }

  iEvent.put(phoCollection);

  if (nSelectedPho >= selection_)
    return true;

  return false;
}


void cicFilter::beginJob()
{}

void cicFilter::endJob() 
{}


float cicFilter::SumPt2(reco::VertexRef vtx) {
    
  float ptsum = 0;

  for(std::vector<reco::TrackBaseRef>::const_iterator trkItr = vtx->tracks_begin();trkItr != vtx->tracks_end(); ++trkItr) {
    ptsum += (*trkItr)->pt();
  }

  return ptsum;
}

int cicFilter::ChooseVertex(edm::Handle<reco::VertexCollection> vtxHandle) {

  int goodIndex = -1;

  float maxPt = 0;

  for(unsigned int i=0; i<vtxHandle->size(); i++) {
    float temp = SumPt2(reco::VertexRef(vtxHandle, i));
    if (temp > maxPt) {
      maxPt = temp;
      goodIndex = i;
    }
  }
  
  return goodIndex;
}


//define this as a plug-in
DEFINE_FWK_MODULE(cicFilter);
