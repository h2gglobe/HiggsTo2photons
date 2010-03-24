#include <string>
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include <iostream>

GlobeJets::GlobeJets(const edm::ParameterSet& iConfig, const char* n = "it5"): nome(n) {
  
  char a[100];
  sprintf (a,"JetColl_%s", nome);
  jetColl =  iConfig.getParameter<edm::InputTag>(a);
  calotowerColl =  iConfig.getParameter<edm::InputTag>("CaloTowerColl");
  trackColl =  iConfig.getParameter<edm::InputTag>("TrackColl");
  
  doEgammaSummer09Skim = iConfig.getParameter<bool>("doEgammaSummer09Skim");
  
  std::string strnome = nome;
  if (strnome.find("pf",0) != std::string::npos){
  }
  
  if (strnome.find("pf",0) == std::string::npos){
    sprintf (a,"JetTrackAssociationColl_%s", nome);
    jetTkAssColl =  iConfig.getParameter<edm::InputTag>(a);
  }

  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeJets::defineBranch(TTree* tree) {

  jet_p4 = new TClonesArray("TLorentzVector", MAX_JETS);
  
  char a1[50], a2[50];
  
  sprintf(a1, "jet_%s_n", nome);
  sprintf(a2, "jet_%s_n/I", nome);
  tree->Branch(a1, &jet_n, a2);
  
  sprintf(a1, "jet_%s_p4", nome);
  tree->Branch(a1, "TClonesArray", &jet_p4, 32000, 0);
  
  sprintf(a1, "jet_%s_emfrac", nome);
  sprintf(a2, "jet_%s_emfrac[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_emfrac, a2);

  sprintf(a1, "jet_%s_hadfrac", nome);
  sprintf(a2, "jet_%s_hadfrac[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_hadfrac, a2);
    
  sprintf(a1, "jet_%s_ntk", nome);
  sprintf(a2, "jet_%s_ntk[jet_%s_n]/I", nome, nome);
  tree->Branch(a1, &jet_ntk, a2);
  
  sprintf(a1, "jet_%s_tkind", nome);
  sprintf(a2, "jet_%s_tkind[jet_%s_n][%d]/I", nome, nome, MAX_JET_TRACKS);
  tree->Branch(a1, &jet_tkind, a2);

  sprintf(a1, "jet_%s_ncalotw", nome);
  sprintf(a2, "jet_%s_ncalotw[jet_%s_n]/I", nome, nome);
  tree->Branch(a1, &jet_ncalotw, a2);
  
  sprintf(a1, "jet_%s_calotwind", nome);
  sprintf(a2, "jet_%s_calotwind[jet_%s_n][%d]/I", nome, nome, MAX_JET_TOWERS);
  tree->Branch(a1, &jet_calotwind, a2);
}

bool GlobeJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::string strnome = nome;
  // if( nome != "pf")
  if (strnome.find("pf",0) == std::string::npos)
    {
      edm::Handle<reco::CaloJetCollection> jetH;
      iEvent.getByLabel(jetColl, jetH);  
      
      // take collections
      
      
      edm::Handle<CaloTowerCollection> ctH; 
      iEvent.getByLabel(calotowerColl, ctH);
      
      edm::Handle<reco::TrackCollection> tkH;
      iEvent.getByLabel(trackColl, tkH);
      
      jet_p4->Clear();
      
      jet_n = 0;
      
      if (debug_level > 9)
	std::cout << "GlobeJets: Jet collection size: "<< jetH->size() << std::endl;
      
      // check if collection is present
      for(unsigned int i=0; i<jetH->size(); i++) {
	if (jet_n >= MAX_JETS) {
	  std::cout << "GlobeJets: WARNING TOO MANY JETS: " << jetH->size() << " (allowed " << MAX_JETS << ")" << std::endl;
	  break;
	}
	
	reco::CaloJetRef j(jetH, i);
	
	// apply the cuts
	if (gCUT->cut(*j)) continue;
	// passed cuts
	
	new ((*jet_p4)[jet_n]) TLorentzVector();
	((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(j->px(), j->py(), j->pz(), j->energy()); 
	jet_emfrac[jet_n] = j->emEnergyFraction();
	jet_hadfrac[jet_n] = j->energyFractionHadronic();
	
	// Tracks and CaloTowers
	if(!doEgammaSummer09Skim) {
	  std::vector<CaloTowerPtr> towers = j->getCaloConstituents();
	  std::vector<CaloTowerPtr>::const_iterator it;
	  
	  int limit = 0;
	  if (towers.size() >= MAX_JET_TOWERS) {
	    std::cout << "GlobeJets: WARNING TOO MANY TOWERS IN JET: " << towers.size() << " (allowed " << MAX_JET_TOWERS << ")" << std::endl;
	    limit = MAX_JET_TOWERS;
	  } else {
	    limit = towers.size();
	  }
	  
	  jet_ncalotw[jet_n] = limit;
	  
	  int index = 0;
	  for(it = towers.begin(); it != towers.end(); ++it) {
	    if (index >= limit)
	      break;
	    for(unsigned int k = 0; k<ctH->size(); k++) {
	      CaloTowerRef t(ctH, k);
	      if (&(**it) == &(*t)) {
		jet_calotwind[jet_n][index] = k;
		break;
	      }
	    }
	    index++;
	  }
	  
	  edm::Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
	  iEvent.getByLabel(jetTkAssColl, jetTracksAssociation);
	  
	  for(reco::JetTracksAssociationCollection::const_iterator itass = jetTracksAssociation->begin(); itass != jetTracksAssociation->end(); ++itass) {
	    if (&(*(itass->first)) != &(*j)) 
	      continue;
	    limit = 0;
	    reco::TrackRefVector tracks = itass->second;
	    if (tracks.size() >= MAX_JET_TRACKS)
	      limit = MAX_JET_TRACKS;
	    else 
	      limit = tracks.size();
	    
	    jet_ntk[jet_n] = limit;
	    
	    for (int ii = 0; ii < limit; ++ii) {
	      for(unsigned int k = 0; k<tkH->size(); k++) {
		reco::TrackRef t(tkH, k);
		if (&(*(tracks[ii])) == &(*t) ) {
		  jet_tkind[jet_n][ii] = k;
		  break;
		}
	      }
	    }
	  }
	}
	jet_n++;
      }
    }
  
  
  if (strnome.find("pf",0) != std::string::npos){
    
    // take collections
    edm::Handle<reco::PFJetCollection> pfjetH;
    iEvent.getByLabel(jetColl, pfjetH);
    
    edm::Handle<CaloTowerCollection> ctH; 
    iEvent.getByLabel(calotowerColl, ctH);
    
    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(trackColl, tkH);
    
    jet_p4->Clear();
    
    jet_n = 0;
    
    if (debug_level > 9)
      std::cout << "GlobeJets: Jet collection size: "<< pfjetH->size() << std::endl;
    
    // check if collection is present
    for(unsigned int i=0; i<pfjetH->size(); i++) {
      if (jet_n >= MAX_JETS) {
	std::cout << "GlobeJets: WARNING TOO MANY JETS: " << pfjetH->size() << " (allowed " << MAX_JETS << ")" << std::endl;
	break;
      }
      
      reco::PFJetRef j(pfjetH, i);
      
      // apply the cuts
      if (gCUT->cut(*j)) continue;
      // passed cuts
      
      new ((*jet_p4)[jet_n]) TLorentzVector();
      ((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(j->px(), j->py(), j->pz(), j->energy()); 
      jet_emfrac[jet_n] = j->chargedEmEnergyFraction() + j->neutralEmEnergyFraction() + j->chargedMuEnergyFraction();
      jet_hadfrac[jet_n] = j->chargedHadronEnergyFraction() + j->neutralHadronEnergyFraction();
      
      jet_n++;
      
    }
    
  }

  return true;
}
