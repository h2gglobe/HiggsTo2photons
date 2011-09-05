// Package:   diphotonFilter 
// Class:     diphotonFilter 
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Jun 22 13:56:58 CEST 2010
// copied/altered version of multiSCfilter.cc, which was a copied/altered version of InvariantMassFilter.cc

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"

class diphotonFilter : public edm::EDFilter {
  public:
    explicit diphotonFilter(const edm::ParameterSet&);
    ~diphotonFilter();

  private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    edm::InputTag phoColl_;

    bool applyFilter_;
    bool applyIso1_;
    bool applyIso2_;
    bool applyID1_;
    bool applyID2_;
    Float_t minEt1_;
    Float_t minEt2_;
    Float_t maxEta1_;
    Float_t maxEta2_;
};


diphotonFilter::diphotonFilter(const edm::ParameterSet& iConfig) {

  phoColl_  = iConfig.getParameter<edm::InputTag>("photonCollection");

  applyFilter_ =  iConfig.getUntrackedParameter("applyFilter",true);
  applyIso1_   =  iConfig.getUntrackedParameter("applyIso1",false);
  applyIso2_   =  iConfig.getUntrackedParameter("applyIso2",false);
  applyID1_    =  iConfig.getUntrackedParameter("applyID1",false);
  applyID2_    =  iConfig.getUntrackedParameter("applyID2",false);
  minEt1_      =  iConfig.getUntrackedParameter("minEt1",0.);
  minEt2_      =  iConfig.getUntrackedParameter("minEt2",0.);
  maxEta1_     =  iConfig.getUntrackedParameter("maxEta1",3.);
  maxEta2_     =  iConfig.getUntrackedParameter("maxEta2",3.);
}

diphotonFilter::~diphotonFilter()
{}

bool diphotonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //bool sel = false;

  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(phoColl_, phoH);

  int nPhotonsPassing=0;



  bool has_good_photon1 = false;
  bool has_good_photon2 = false;

  if(!applyFilter_)return true;
  //if(phoH->size() < 2)std::cout << "too few photons in event" << std::endl;
  if(phoH->size() < 2)return false;

  for(reco::PhotonCollection::const_iterator pho = phoH->begin(); pho != phoH->end(); pho++) {

    bool pass_ecal_iso = pho->ecalRecHitSumEtConeDR04() < 1.5*(4.2+0.006*pho->pt());
    bool pass_hcal_iso = pho->hcalTowerSumEtConeDR04()  < 1.5*(2.2+0.0025*pho->pt());
    bool pass_trk_iso  = pho->trkSumPtHollowConeDR04()  < 1.5*(3.5+0.001*pho->pt());
    bool pass_sieie_id = pho->sigmaIetaIeta()  < 0.013*pho->isEB() + 0.03*pho->isEE();
    bool pass_hoe_id   = pho->hadronicOverEm()  < 0.05;
    bool pass_pixel_veto   = true;//!(pho->hasPixelSeed());

    bool thisPhotonPassesIsolation1 = !applyIso1_ || (pass_ecal_iso && pass_hcal_iso && pass_trk_iso);
    bool thisPhotonPassesID1 = !applyID1_ || (pass_sieie_id && pass_hoe_id && pass_pixel_veto);
    bool thisPhotonPassesEtaCut1 = (fabs(pho->caloPosition().eta()) < maxEta1_);
    bool thisPhotonPassesEtCut1 = (sin(pho->caloPosition().theta())*pho->superCluster()->energy() > minEt1_);
    bool thisPhotonPassedCuts1 =  thisPhotonPassesIsolation1 && thisPhotonPassesID1 && thisPhotonPassesEtaCut1 && thisPhotonPassesEtCut1;
    //std::cout << "photon1 iso,id,eta,et: " << thisPhotonPassesIsolation1 <<","<< thisPhotonPassesID1 <<","<< thisPhotonPassesEtaCut1 <<","<< thisPhotonPassesEtCut1 << std::endl;

    bool thisPhotonPassesIsolation2 = !applyIso2_ || (pass_ecal_iso && pass_hcal_iso && pass_trk_iso);
    bool thisPhotonPassesID2 = !applyID2_ || (pass_sieie_id && pass_hoe_id && pass_pixel_veto);
    bool thisPhotonPassesEtaCut2 = (fabs(pho->caloPosition().eta()) < maxEta2_);
    bool thisPhotonPassesEtCut2 = (sin(pho->caloPosition().theta())*pho->superCluster()->energy() > minEt2_);
    bool thisPhotonPassedCuts2 = thisPhotonPassesIsolation2 && thisPhotonPassesID2 && thisPhotonPassesEtaCut2 && thisPhotonPassesEtCut2;
    //std::cout << "photon2 iso,id,eta,et: " << thisPhotonPassesIsolation2 <<","<< thisPhotonPassesID2 <<","<< thisPhotonPassesEtaCut2 <<","<< thisPhotonPassesEtCut2 << std::endl;


    // keep track if this event has a photon passing photon cuts 1
    if(thisPhotonPassedCuts1) has_good_photon1=true;
    // keep track if this event has a photon passing photon cuts 2
    if(thisPhotonPassedCuts2) has_good_photon2=true;
    // count how many photons have passed either cut 
    if(thisPhotonPassedCuts1 || thisPhotonPassedCuts2) {
      nPhotonsPassing++;
    } else {
      continue;
    }

    //if(thisPhotonPassedCuts1)std::cout << "PASSED PHOTON 1" << std::endl;
    //if(thisPhotonPassedCuts2)std::cout << "PASSED PHOTON 2" << std::endl;

    // keep event if have found a photon passing cuts 1, a photon passing cuts 2, and at least TWO photons total.
    if(has_good_photon1 && has_good_photon2 && nPhotonsPassing > 1)return true;

  }
  return false;
}

void diphotonFilter::beginJob()
{}

void diphotonFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(diphotonFilter);
