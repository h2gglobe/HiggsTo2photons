#ifndef GLOBECUTS_H
#define GLOBECUTS_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsTo2photons//interface/Limits.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Math/VectorUtil.h"

#include <iostream>

class GlobeCuts {
 public:
  
  GlobeCuts(const edm::ParameterSet&);
  virtual ~GlobeCuts() {};

  bool cut(const reco::Photon&);
  bool cut(const reco::PFCandidate&);
  bool cut(const reco::GsfElectron&);
  bool cut(const reco::SuperCluster&);
  bool cut(const reco::BasicCluster&);
  bool cut(const CaloTower&);
  bool cut(const HBHERecHit&,double);
  bool cut(const HFRecHit&);
  bool cut(const HORecHit&);
  bool cut(const reco::Track&);
  bool cut(const PSimHit&);
  bool cut(const SimTrack&);
  bool cut(const reco::CaloJet&);
  bool cut(const reco::PFJet&);
  bool cut(const reco::BasicJet&);
  bool cut(const TrackingParticle&);
  bool cut(const reco::Muon&);
  bool cut(const reco::GenJet&);
  bool cut(const EcalRecHit&,int,double);

  bool isocut(const reco::Track&,const reco::Track&,const reco::Vertex&);
  bool isocut(const reco::Track&,const reco::Track&);

 private:
 /*
  edm::ParameterSet psetPhoton;
  edm::ParameterSet psetConvertedPhoton;
  edm::ParameterSet psetElectron;
  edm::ParameterSet psetSuperCluster;
  edm::ParameterSet psetBasicCluster;
  edm::ParameterSet psetCaloTower;
  edm::ParameterSet psetHcalHit;
  edm::ParameterSet psetTrack;
  edm::ParameterSet psetSimTrack;
  edm::ParameterSet psetJet;
  edm::ParameterSet psetMuon;
  edm::ParameterSet psetGenJet;
  edm::ParameterSet psetEcalHit;
  edm::ParameterSet psetIso;
  edm::ParameterSet psetTP;
*/

  double photonEtCut_; 
  double convertedPhotonEtCut_; 
  double electronEtCut_; 
  double superClusterEnergyCut_; 
  double basicClusterEnergyCut_; 
  double caloTowerEtCut_; 

  double hcalHBHEEnergyCut_; 
  double hcalMaxDR_; 
  double hcalHFEnergyCut_; 
  double hcalHOEnergyCut_; 
  bool   hcalKeepOutsideCone_;

  std::vector<int> tpPdgidCut_;
  double trackPtCut_; 
  double tpPtCut_;
  double tpEtaCut_;
  double tpTvpCut_;
  double tpLvpCut_;
  double simHitEnergyCut_; 
  double simTrackEnergyCut_; 
  double jetEnergyCut_; 
  double muonPtCut_; 
  double genJetEtCut_;
 
  double ecalHitBarrelECut_;
  double ecalHitEndcapECut_;
  double ecalHitPreECut_;
  bool   ecalKeepOutsideCone_;
  double ecalMaxDR_;

  double isoInnerCone_;
  double isoOuterCone_;

};


#endif
