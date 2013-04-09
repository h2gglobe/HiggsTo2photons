#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"

#include <iostream>

GlobeGenParticles::GlobeGenParticles(const edm::ParameterSet& iConfig) {
  
  genParticlesColl = iConfig.getParameter<edm::InputTag>("GenParticlesColl");
  debug_level = iConfig.getParameter<int>("Debug_Level");
  gCUT = new GlobeCuts(iConfig);
}

void GlobeGenParticles::defineBranch(TTree* tree) {

  // think about changing branch names for duplicate collections
  gp_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
  gp_vtx = new TClonesArray("TVector3", MAX_GENERATOR);
  
  tree->Branch("gp_n", &gp_n, "gp_n/I");
  
  tree->Branch("gp_p4", "TClonesArray", &gp_p4, 32000, 0);
  tree->Branch("gp_vtx", "TClonesArray", &gp_vtx, 32000, 0);
  
  tree->Branch("gp_status", gp_status, "gp_status[gp_n]/S");
  tree->Branch("gp_pdgid", gp_pdgid, "gp_pdgid[gp_n]/S");
  tree->Branch("gp_mother", gp_mother, "gp_mother[gp_n]/S");
}

bool GlobeGenParticles::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GenParticleCollection> gpH;
  iEvent.getByLabel(genParticlesColl, gpH);

  gp_p4->Clear();
  gp_vtx->Clear();
  gp_n = 0;
  
  if(debug_level > 999)
    std::cout << "#"
	      << "\t pdgid"
	      << "\t status"
	      << "\t mother"
	      << "\t px \t py \t pz \t E"
	      << "\t vx \t vy \t vz"
	      << std::endl;
  
  
  for(size_t i = 0; i < gpH->size(); ++i) {

    const reco::GenParticleRef gp(gpH, i);

    if(!gCUT->cut(*gp))
      continue;

    gp_pdgid[gp_n] = gp->pdgId();
    gp_status[gp_n] = gp->status();
    if (gp->numberOfMothers() != 0)
      gp_mother[gp_n] = gp->motherRef().key();
    else
      gp_mother[gp_n] = -1;

    new ((*gp_p4)[gp_n]) TLorentzVector();
    ((TLorentzVector *)gp_p4->At(gp_n))->SetPtEtaPhiM(gp->pt(), gp->eta(), gp->phi(), gp->mass());

    new ((*gp_vtx)[gp_n]) TVector3();
    ((TVector3 *)gp_vtx->At(gp_n))->SetXYZ(gp->vx(), gp->vy(), gp->vz());
    
    if(debug_level > 999)
      std::cout << gp_n	<< "\t" << gp_pdgid[ gp_n]
                << "\t" << gp_status[gp_n] << "\t" << gp_mother[gp_n]
                << "\t" << gp->px() << "\t" << gp->py() 
                << "\t" << gp->pz() << "\t" << gp->energy()
                << "\t" << gp->vx() << "\t" << gp->vy() << "\t" << gp->vz()
                << std::endl;
    
    gp_n++;
  }
  
  return true;
}
