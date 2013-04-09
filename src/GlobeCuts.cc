#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

GlobeCuts::GlobeCuts(const edm::ParameterSet& iConfig) {

  edm::ParameterSet psetGenParticle         = iConfig.getParameter<edm::ParameterSet>("GenParticleCuts");
  edm::ParameterSet psetPhoton              = iConfig.getParameter<edm::ParameterSet>("PhotonCuts");
  edm::ParameterSet psetConvertedPhoton     = iConfig.getParameter<edm::ParameterSet>("ConvertedPhotonCuts");
  edm::ParameterSet psetElectron            = iConfig.getParameter<edm::ParameterSet>("ElectronCuts");
  edm::ParameterSet psetSuperCluster        = iConfig.getParameter<edm::ParameterSet>("SuperClusterCuts");
  edm::ParameterSet psetBasicCluster        = iConfig.getParameter<edm::ParameterSet>("BasicClusterCuts");
  edm::ParameterSet psetCaloTower           = iConfig.getParameter<edm::ParameterSet>("CaloTowerCuts");
  edm::ParameterSet psetHcalHit             = iConfig.getParameter<edm::ParameterSet>("HcalHitsCuts");
  edm::ParameterSet psetTrack               = iConfig.getParameter<edm::ParameterSet>("TrackCuts");
  edm::ParameterSet psetSimHit              = iConfig.getParameter<edm::ParameterSet>("SimHitCuts");
  edm::ParameterSet psetSimTrack            = iConfig.getParameter<edm::ParameterSet>("SimTrackCuts");
  edm::ParameterSet psetJet                 = iConfig.getParameter<edm::ParameterSet>("JetCuts");
  edm::ParameterSet psetMuon                = iConfig.getParameter<edm::ParameterSet>("MuonCuts");
  edm::ParameterSet psetGenJet              = iConfig.getParameter<edm::ParameterSet>("GenJetCuts");
  edm::ParameterSet psetEcalHit             = iConfig.getParameter<edm::ParameterSet>("EcalHitCuts");
  edm::ParameterSet psetIso                 = iConfig.getParameter<edm::ParameterSet>("IsoCuts");
  edm::ParameterSet psetTP                  = iConfig.getParameter<edm::ParameterSet>("TPCuts");

  genParticlePdgidCut_      = psetGenParticle.getParameter<std::vector<int> >("PdgId");
  genParticleBoolCut_       = psetGenParticle.getParameter<bool>("Keep");
  genParticleEtCut_         = psetGenParticle.getParameter<double>("EtCut");
  photonEtCut_              = psetPhoton.getParameter<double>("EtCut");
  convertedPhotonEtCut_     = psetConvertedPhoton.getParameter<double>("EtCut");
  electronEtCut_            = psetElectron.getParameter<double>("EtCut");
  superClusterEnergyCut_    = psetSuperCluster.getParameter<double>("EnergyCut");
  basicClusterEnergyCut_    = psetBasicCluster.getParameter<double>("EnergyCut");
  caloTowerEtCut_           = psetCaloTower.getParameter<double>("EtCut");
  trackPtCut_               = psetTrack.getParameter<double>("PtCut");
  simHitEnergyCut_          = psetSimHit.getParameter<double>("EnergyCut");
  simTrackEnergyCut_        = psetSimTrack.getParameter<double>("EnergyCut");
  jetEnergyCut_             = psetJet.getParameter<double>("EnergyCut");
  muonPtCut_                = psetMuon.getParameter<double>("PtCut");
  genJetEtCut_              = psetGenJet.getParameter<double>("EtCut");

  tpPtCut_                  = psetTP.getParameter<double>("tpPtCut");
  tpEtaCut_                 = psetTP.getParameter<double>("tpEtaCut");
  tpTvpCut_                 = psetTP.getParameter<double>("tpTvpCut");
  tpLvpCut_                 = psetTP.getParameter<double>("tpLvpCut");
  tpPdgidCut_               = psetTP.getParameter<std::vector<int> > ("tpPdgidCut");

  hcalHBHEEnergyCut_        = psetHcalHit.getParameter<double>("HBHEEnergyCut");
  hcalMaxDR_                = psetHcalHit.getParameter<double>("HcalMaxDR");
  hcalHFEnergyCut_	    = psetHcalHit.getParameter<double>("HFEnergyCut");
  hcalKeepOutsideCone_      = psetEcalHit.getParameter<bool>("KeepOutsideCone");
  hcalHOEnergyCut_	    = psetHcalHit.getParameter<double>("HOEnergyCut");

  ecalHitBarrelECut_        = psetEcalHit.getParameter<double>("BarrelEnergyCut");
  ecalHitEndcapECut_        = psetEcalHit.getParameter<double>("EndcapEnergyCut");
  ecalHitPreECut_           = psetEcalHit.getParameter<double>("PreEnergyCut");
  ecalKeepOutsideCone_      = psetEcalHit.getParameter<bool>("KeepOutsideCone");
  ecalMaxDR_                = psetEcalHit.getParameter<double>("EcalMaxDR");

  isoInnerCone_ 	    = psetIso.getParameter<double>("InnerCone");
  isoOuterCone_ 	    = psetIso.getParameter<double>("OuterCone");
}

// The Functions return "true" if the object should be cut
// GenParticle
bool GlobeCuts::cut(const reco::GenParticle &gp) {

  if (gp.et() < genParticleEtCut_)
    return 0;

  if (genParticlePdgidCut_.size() == 0)
    return 1;

  for(unsigned int i=0; i<genParticlePdgidCut_.size(); i++)
    if ((abs(gp.pdgId()) == genParticlePdgidCut_[i] and genParticleBoolCut_) or 
	(abs(gp.pdgId()) != genParticlePdgidCut_[i] and !genParticleBoolCut_))
      return 1;

  return 0;
}

//PFCands
bool GlobeCuts::cut(const reco::PFCandidate &pf) {
  return 0;
}

// Photons
bool GlobeCuts::cut(const reco::Photon &photon) { 
  return (photon.et() < photonEtCut_); 
}


// Conversions
bool GlobeCuts::cut(const reco::Conversion &conv) { 
  return (sqrt(conv.refittedPairMomentum().perp2()) <  convertedPhotonEtCut_); 
}


// Electrons
bool GlobeCuts::cut(const reco::GsfElectron& electron) { 
          return (electron.et() < electronEtCut_); 
}
// Super Clusters
bool GlobeCuts::cut(const reco::SuperCluster &supercluster) { 
          return (supercluster.energy() < superClusterEnergyCut_); 
}
// Basic Clusters
bool GlobeCuts::cut(const reco::BasicCluster &basiccluster) { 
          return (basiccluster.energy() < basicClusterEnergyCut_); 
}
// CaloTowers
bool GlobeCuts::cut(const CaloTower &calotower) { 
          return (calotower.et() < caloTowerEtCut_); 
}

// Hcal RecHits Barrel
bool GlobeCuts::cut(const HBHERecHit &hcalrechit, double dR) { 

    if( dR < hcalMaxDR_ ) return false; //keep if inside cone

    if (!hcalKeepOutsideCone_) return true; //dont keep if user doesnt want

    return (hcalrechit.energy() < hcalHBHEEnergyCut_ );
}

// Hcal RecHits Forward
bool GlobeCuts::cut(const HFRecHit &hcalrechit) { 
          return (hcalrechit.energy() < hcalHFEnergyCut_); 
}

// Hcal RecHits Outer
bool GlobeCuts::cut(const HORecHit &hcalrechit) { 
          return (hcalrechit.energy() < hcalHOEnergyCut_); 
}

// Tracks
bool GlobeCuts::cut(const reco::Track &track) { 
          return (track.pt() < trackPtCut_); 
}

// TrackingParticle
bool GlobeCuts::cut(const TrackingParticle &tp) { 

  bool pdgidCut = true;
  for(unsigned int i=0; i<tpPdgidCut_.size(); i++) {
    if (abs(tp.pdgId()) == tpPdgidCut_[i]) {
      pdgidCut = false;
      break;
    }
  }

  return ( tp.pt() < tpPtCut_ || 
           fabs(tp.eta()) > tpEtaCut_ || 
           fabs(tp.vertex().z()) > tpLvpCut_ || 
           fabs(tp.vertex().rho()) > tpTvpCut_ ||
           pdgidCut); 
}


// Sim Hits
bool GlobeCuts::cut(const PSimHit &simhit) { 
          return (simhit.energyLoss() < simHitEnergyCut_); 
}
// Sim Tracks
bool GlobeCuts::cut(const SimTrack &simtrack) { 
          return (simtrack.momentum().E() < simTrackEnergyCut_ || simtrack.noVertex());
}
// Jets
bool GlobeCuts::cut(const reco::CaloJet &jet) { 
  return (jet.energy() < jetEnergyCut_); 
}

bool GlobeCuts::cut(const reco::PFJet &jet) { 
  return (jet.energy() < jetEnergyCut_); 
}

bool GlobeCuts::cut(const reco::BasicJet &jet) { 
  return (jet.energy() < jetEnergyCut_); 
}

// Muons
bool GlobeCuts::cut(const reco::Muon &muon) { 
          return (muon.pt() < muonPtCut_); 
}
// Gen Jets
bool GlobeCuts::cut(const reco::GenJet &genjet) { 
          return (genjet.et() < genJetEtCut_); 
}

//Cut to determine if RecHit falls within the lepton cone
//const EcalRecHit here is just used as a function identifier
bool GlobeCuts::cut(const EcalRecHit &ecalhit, int type, double dR) { 

  if (type == 2)  // preshower
    return (ecalhit.energy() < ecalHitPreECut_); 
   
  if( dR < ecalMaxDR_ ) 
    return false; //keep if inside cone
  
  if (!ecalKeepOutsideCone_) 
    return true; //dont keep if user doesnt want stuff outside cone

  if (type == 0)  // barrel
    return (fabs(ecalhit.energy()) < ecalHitBarrelECut_); 
  if (type == 1)  // endcap
    return (fabs(ecalhit.energy()) < ecalHitEndcapECut_); 
  return false;
}

bool GlobeCuts::isocut(const reco::Track &tk, const reco::Track &lep, const reco::Vertex &vtx) {
    
    using ROOT::Math::VectorUtil::DeltaR;
    
    //Inside the veto cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) < isoInnerCone_) return false;

    //Outside the isolation cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) > isoOuterCone_) return false;

    double d0 = sqrt((tk.vertex() - vtx.position()).perp2());
    double dZ = (tk.vertex() - vtx.position()).z();

    if(tk.numberOfValidHits() >= 10)
        return (tk.pt() > 1 && d0 < 1 && dZ < 5);
    else if(tk.numberOfValidHits() >= 8 && tk.numberOfValidHits() <= 9)
        return (tk.pt() > 1 && d0 < 0.2 && dZ < 2.0 && (d0 / tk.d0Error() < 10) && (dZ / tk.dzError() < 10));
    else if(tk.numberOfValidHits() >= 5 && tk.numberOfValidHits() <= 7)
        return (tk.pt() > 1 && d0 < 0.04 && dZ < 0.5 && (d0 / tk.d0Error() < 7) && (dZ / tk.dzError() < 7));
    else if(tk.numberOfValidHits() >= 3 && tk.numberOfValidHits() <= 4)
        return (tk.pt() > 1 && d0 < 0.02 && dZ < 0.2 && (d0 / tk.d0Error() < 3) && (dZ / tk.dzError() < 3));
    else return false;
    

}

bool GlobeCuts::isocut(const reco::Track &tk, const reco::Track &lep) {
    
    using ROOT::Math::VectorUtil::DeltaR;

    //Inside the veto cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) < isoInnerCone_) return false;

    //Outside the isolation cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) > isoOuterCone_) return false;

    double d0 = sqrt((tk.vertex() - lep.vertex()).perp2());
    double dZ = (tk.vertex() - lep.vertex()).z();
    
    if(tk.numberOfValidHits() >= 10)
        return (tk.pt() > 1 && d0 < 1 && dZ < 5);
    else if(tk.numberOfValidHits() >= 8 && tk.numberOfValidHits() <= 9)
        return (tk.pt() > 1 && d0 < 0.2 && dZ < 2.0 && d0 / tk.d0Error() < 10 && dZ / tk.dzError() < 10);
    else if(tk.numberOfValidHits() >= 5 && tk.numberOfValidHits() <= 7)
        return (tk.pt() > 1 && d0 < 0.04 && dZ < 0.5 && d0 / tk.d0Error() < 7 && dZ / tk.dzError() < 7);
    else if(tk.numberOfValidHits() >= 3 && tk.numberOfValidHits() <= 4)
        return (tk.pt() > 1 && d0 < 0.02 && dZ < 0.2 && d0 / tk.d0Error() < 3 && dZ / tk.dzError() < 3);
    else return false;
    

}

    




