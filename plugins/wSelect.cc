// Package:    WSelect
// Class:      WSelect
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

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/eIDCuts.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class WSelect : public edm::EDFilter {
public:
  explicit WSelect(const edm::ParameterSet&);
  ~WSelect();

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
  edm::InputTag eIDLabel;
};


WSelect::WSelect(const edm::ParameterSet& iConfig) {
  
  jetColl_ = iConfig.getParameter<edm::InputTag>("jetCollection");
  metColl_ = iConfig.getParameter<edm::InputTag>("metCollection");
  elColl_  = iConfig.getParameter<edm::InputTag>("electronCollection");
  scColl_  = iConfig.getParameter<edm::InputTag>("superClusterCollection");
  
  useSuperCluster_ = iConfig.getParameter<bool>("useSuperCluster");
  metCut_ =  iConfig.getParameter<double>("metCut");
  
  eIDLabel = iConfig.getParameter<edm::InputTag>("eIDLabel");
}

WSelect::~WSelect()
{}

bool WSelect::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::CaloJetCollection> jetH;
  iEvent.getByLabel(jetColl_, jetH);

  //edm::Handle<reco::PFJetCollection> jetH;
  //iEvent.getByLabel(jetColl_, jetH);

  edm::Handle<reco::PFMETCollection> metH;
  iEvent.getByLabel(metColl_, metH);
    
  edm::Handle<edm::ValueMap<float> > eIDVM; 
  iEvent.getByLabel(eIDLabel, eIDVM);
  
//  const edm::ValueMap<float>& eIDmapTemp = *eIDVM;
  Int_t eIso_medium = -1; 

  Float_t etCand = 0;
  math::XYZVector dir(0,0,0);

  if (useSuperCluster_) {
    edm::Handle<reco::SuperClusterCollection> scH;
    iEvent.getByLabel(scColl_, scH);
    dir = scH->begin()->position();
    etCand = scH->begin()->energy() * sin(scH->begin()->position().theta());

  } else {
    edm::Handle<reco::GsfElectronCollection> elH;
    iEvent.getByLabel(elColl_, elH);
    dir =  math::XYZVector(elH->begin()->px(), 
                           elH->begin()->py(), 
                           elH->begin()->pz()); 
    etCand = elH->begin()->et();
    reco::GsfElectronRef electronRef(elH, 0);
    //reco::GsfElectronRef electronRef(elH, std::distance(elH->begin(), elH->begin()));
    //eID_medium = (int)eIDmapTemp[electronRef];
    reco::GsfElectron egsf = reco::GsfElectron(*electronRef);
    
    // Variables for Catorgories
    Float_t fbrem = egsf.fbrem();
    Float_t eopin = egsf.eSuperClusterOverP();

    Int_t cat = -1;
    Float_t eta = egsf.superCluster()->eta();
    if (eta < 0) eta = -eta;

    // Catorgories
    if (eta < 1.479) {
      if (fbrem>=0.12 && eopin >0.9 && eopin < 1.2)    // bremming barrel
        cat = 0;
      else if (fbrem < 0.12)                           // low brem barrel 
        cat = 1;
      else                                             // bad track in barrel
        cat = 2;
    } else {
      if (fbrem>=0.2 && eopin >0.82 && eopin < 1.22)   // bremming endcap
        cat = 3;
      else if (eta > 1.5 && eta < 1.58)                // endcap crack electrons
        cat = 7;                                       
      else if (fbrem < 0.2)                            // low brem endcap 
        cat = 4;
      else                                             // bad track in barrel
        cat = 5;
    }
    
    Int_t bin = 0;

    //  ET bin
    if (etCand < 20.) 
      bin = 2;
    else if (etCand > 30.)
      bin = 0;
    else
      bin = 1;

    
    // Variables to cut
    Float_t el_tkiso03 = egsf.dr03TkSumPt();
    Float_t el_ecaliso04 = egsf.dr04EcalRecHitSumEt();
    Float_t el_hcaliso04 = egsf.dr04HcalTowerSumEt();
    
    //  isosum
    Float_t iso_sum = el_tkiso03 + el_ecaliso04 + el_hcaliso04;   
    iso_sum-=0.15+0.1*(2);                                  //  ad hoc correction of data, nvert = 2  isosum
    //iso_sum-=0.15+0.1*(nvert-1);                          //  ad hoc correction of data             isosum
    if(fabs(eta)>1.479) iso_sum+=0.02;                  //  ad hoc correction of data
    Float_t iso_sum_scaled = iso_sum*pow(40/etCand, 2); 


    Int_t type = 2;  // eIDlevel is medium
    if ((iso_sum < cutiso_sum[bin][type][cat]) &&                                          // isosum
        (iso_sum_scaled < cutiso_sumoet[bin][type][cat]))                                  // sumoet
      eIso_medium = 1;
    else  eIso_medium = 0;

  }

  // MET
  Float_t met = metH->begin()->et();

  bool Wsel = etCand > 18 and met > 18 and eIso_medium == 1; //checking for medium isolation
  //bool Wsel = etCand > 18 and met > 18 and (eID_medium == 15);  

    
  return Wsel;
}

void WSelect::beginJob()
{}

void WSelect::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(WSelect);
