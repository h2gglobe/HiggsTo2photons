#ifndef __PFIsolation__
#define __PFIsolation__

#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollection.h"
#include <vector>

bool vetoPFParticle(const reco::PFCandidate&, std::vector<reco::PFCandidate::ParticleType>);

float pfHcalIso(const reco::GsfElectron&, const reco::PFCandidateCollection*, float, float, std::vector<reco::PFCandidate::ParticleType>);

float pfEcalIso(const reco::GsfElectron&, const reco::PFCandidateCollection*, float, float, std::vector<reco::PFCandidate::ParticleType>);

float pfTkIso(const reco::GsfElectron&, edm::Handle<reco::PFCandidateCollection>, edm::Handle<reco::PileUpPFCandidateCollection>, float, float, std::vector<reco::PFCandidate::ParticleType>);

float hoeCalculator(const reco::BasicCluster*, const CaloGeometry&, edm::Handle<HBHERecHitCollection>);
#endif
