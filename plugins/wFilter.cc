// Package:    WFilter
// Class:      WFilter
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

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class WFilter : public edm::EDFilter {
public:
  explicit WFilter(const edm::ParameterSet&);
  ~WFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag jetColl_;
  edm::InputTag metColl_;
  edm::InputTag elColl_;
  edm::InputTag scColl_;
  
  bool useSuperCluster_;
  Float_t metCut_;
  edm::InputTag eIDLabelCiC, eIDLabelVBTF;
};


WFilter::WFilter(const edm::ParameterSet& iConfig) {
  
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetCollection");
  metColl_ = iConfig.getParameter<edm::InputTag>("metCollection");
  elColl_  = iConfig.getParameter<edm::InputTag>("electronCollection");
  scColl_  = iConfig.getParameter<edm::InputTag>("superClusterCollection");
  
  useSuperCluster_ = iConfig.getParameter<bool>("useSuperCluster");
  metCut_ =  iConfig.getParameter<double>("metCut");
  
  eIDLabelCiC = iConfig.getParameter<edm::InputTag>("eIDLabelCiC");
  eIDLabelVBTF = iConfig.getParameter<edm::InputTag>("eIDLabelVBTF");
}

WFilter::~WFilter()
{}

bool WFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  bool sel = false;
  
  //edm::Handle<reco::CaloJetCollection> jetH;
  //iEvent.getByLabel(jetColl_, jetH);

  edm::Handle<reco::PFJetCollection> jetH;
  iEvent.getByLabel(jetColl_, jetH);

  edm::Handle<reco::PFMETCollection> metH;
  iEvent.getByLabel(metColl_, metH);
    
  Int_t eID_medium = -1; 
  Int_t eID_vbtf80 = -1; 
  bool selEID = false;
 
  Float_t etCand = 0;
  math::XYZVector dir(0,0,0);

  if (useSuperCluster_) {
    edm::Handle<reco::SuperClusterCollection> scH;
    iEvent.getByLabel(scColl_, scH);
    dir = scH->begin()->position();
    etCand = scH->begin()->energy() * sin(scH->begin()->position().theta());
    selEID = true;

  } else {

    edm::Handle<edm::ValueMap<float> > eIDCiCVM; 
    edm::Handle<edm::ValueMap<float> > eIDVBTFVM; 
    
    edm::Handle<reco::GsfElectronCollection> elH;
    iEvent.getByLabel(elColl_, elH);
    dir =  math::XYZVector(elH->begin()->px(), 
                           elH->begin()->py(), 
                           elH->begin()->pz()); 
    etCand = elH->begin()->et();
    reco::GsfElectronRef electronRef(elH, 0);

    
    if (eIDLabelCiC.label() != "") { 
      iEvent.getByLabel(eIDLabelCiC, eIDCiCVM);
      const edm::ValueMap<float>& eIDmapTempCiC = *eIDCiCVM;
      eID_medium = (int)eIDmapTempCiC[electronRef];
    } else {
      eID_medium = 15;
    }
    
    if (eIDLabelVBTF.label() != "") {
      iEvent.getByLabel(eIDLabelVBTF, eIDVBTFVM);
      const edm::ValueMap<float>& eIDmapTempVBTF = *eIDVBTFVM; 
      eID_vbtf80 = (int)eIDmapTempVBTF[electronRef];
    } else {
      eID_vbtf80 = 7;
    }

    if (eID_medium == 15 || eID_vbtf80 == 7)
      selEID = true;
  }

  // MET
  Float_t met = metH->begin()->et();

  // dPhi candidate MET
  Float_t dPhi = fabs(deltaPhi(metH->begin()->phi(), dir.phi()));
  
  // sumET
  Float_t sumEt = 0;
  for(unsigned int i=0; i<jetH->size(); i++) {
    
    //reco::CaloJetRef j(jetH, i);
    reco::PFJetRef j(jetH, i);
    if (j->energy() > 5.) {
      Float_t dR = deltaR(j->p4(), dir);
      if (dR > 0.4)
        sumEt += j->et();
    }
  }
  
  bool sel0 = etCand > 18 and etCand < 60 and dPhi > 0.75;
  // new cut using pf jets
  bool sel30 = (met > (metCut_ - 10*(met/sumEt-2))) and ((met/sumEt) > 0.3) and sel0;  
  // old cut using calo jets
  //bool sel30 = (met > (metCut_ - 10*(met/sumEt-2))) and ((met/sumEt) > 0.5) and sel0;  
  sel = sel30 and !(sumEt<10 and met<(20+0.8*sumEt)) and selEID;
    
  return sel;
}

void WFilter::beginJob()
{}

void WFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(WFilter);
