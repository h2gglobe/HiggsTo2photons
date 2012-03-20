#include "HiggsAnalysis/HiggsTo2photons/interface/PFIsolation.h"

#include "DataFormats/Math/interface/deltaR.h"

bool vetoPFParticle(const reco::PFCandidate& pfc, std::vector<reco::PFCandidate::ParticleType> pVetoes) {

  for (std::vector<reco::PFCandidate::ParticleType>::const_iterator it = pVetoes.begin(); it != pVetoes.end(); ++it) {
    if (pfc.particleId() == *it) {
      return false;
    }
  }

  return true;
}    


float pfHcalIso(const reco::GsfElectron& egsf, const reco::PFCandidateCollection* pfCollection, float dRmax, float dRveto, std::vector<reco::PFCandidate::ParticleType> pVetoes) {

  math::XYZVector dir(egsf.px(), egsf.py(), egsf.pz());

  float sum = 0;
  for(unsigned i=0; i<pfCollection->size(); i++) {
    
    const reco::PFCandidate& pfc = (*pfCollection)[i];

    if (!vetoPFParticle(pfc, pVetoes)) {

      reco::GsfTrackRef pfGsfTrackRef = pfc.gsfTrackRef();
      if(egsf.gsfTrack() == pfGsfTrackRef) 
	continue;
	
      math::XYZVector pvi(pfc.momentum());
      float dR = deltaR(dir.Eta(), dir.Phi(), pvi.Eta(), pvi.Phi());
      
      if(dR > dRmax || dR < dRveto)
	continue;

      sum += pfc.pt();
    }
  }
  
  return sum;
}

float pfEcalIso(const reco::GsfElectron& egsf, const reco::PFCandidateCollection* pfCollection, float dRmax, float dRveto, std::vector<reco::PFCandidate::ParticleType> pVetoes) {

  math::XYZVector dir(egsf.px(), egsf.py(), egsf.pz());
  float sum = 0;
  for(unsigned i=0; i<pfCollection->size(); i++) {
    
    const reco::PFCandidate& pfc = (*pfCollection)[i];

    if (!vetoPFParticle(pfc, pVetoes)) {
      math::XYZVector pvi(pfc.momentum());
      float dR = deltaR(dir.Eta(), dir.Phi(), pvi.Eta(), pvi.Phi());
      
      if(dR > dRmax || dR < dRveto)
	continue;
      
      sum += pfc.pt();
    }
  }
  
  return sum;
}

float pfTkIso(const reco::GsfElectron& egsf, edm::Handle<reco::PFCandidateCollection> pfCollectionH, edm::Handle<reco::PFCandidateCollection> pfCollectionPUH, float dRmax, float dRveto, std::vector<reco::PFCandidate::ParticleType> pVetoes) {

  math::XYZVector dir(egsf.px(), egsf.py(), egsf.pz());
  edm::ProductID originalID = pfCollectionH.id();
  float sum = 0;

  for(unsigned i=0; i<pfCollectionH->size(); i++) {
    edm::Ptr<reco::PFCandidate> pfc(pfCollectionH, i);
    
    if (!vetoPFParticle(*pfc, pVetoes)) {

      bool isCandFromPU = false;
      for (unsigned npuit=0; npuit < pfCollectionPUH->size(); npuit++)   {
	
	edm::Ptr<reco::PFCandidate>  ptr(pfCollectionPUH, npuit);
	reco::CandidatePtr candPtr(ptr);
	unsigned nSources = candPtr->numberOfSourceCandidatePtrs();
	if(nSources > 1) {
	  for(unsigned i=0; i<nSources; i++) {
	    reco::CandidatePtr mother = candPtr->sourceCandidatePtr(i);
	    if(originalID == mother.id()) {
	      if(mother.key() == i) {
		isCandFromPU = true;
		break;
	      }
	    }
	  }
	}
      }
      
      if (!isCandFromPU) {
	math::XYZVector pvi(pfc->momentum());
	float dR = deltaR(dir.Eta(), dir.Phi(), pvi.Eta(), pvi.Phi());
	
	if(dR > dRmax || dR < dRveto)
	  continue;
	
	sum += pfc->pt();
      }
    }
  }

  return sum;
}

