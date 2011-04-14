#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenVertices.h"

#include <iostream>

using namespace edm;

GlobeGenVertices::GlobeGenVertices(const edm::ParameterSet& iConfig) {
  
  genParticlesColl = iConfig.getParameter<edm::InputTag>("GenParticlesColl");
  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeGenVertices::defineBranch(TTree* tree) {

  // think about changing branch names for duplicate collections
  gv_pos = new TClonesArray("TVector3", MAX_VERTICES);
  gv_p3 = new TClonesArray("TVector3", MAX_VERTICES);
  
  tree->Branch("gv_n", &gv_n, "gv_n/I");
  
  tree->Branch("gv_pos", "TClonesArray", &gv_pos, 32000, 0);
  tree->Branch("gv_p3", "TClonesArray", &gv_p3, 32000, 0);
  
  tree->Branch("gv_sumPtHi", gv_sumPtHi, "gv_sumPtHi[gv_n]/F");
  tree->Branch("gv_sumPtLo", gv_sumPtLo, "gv_sumPtLo[gv_n]/F");
  tree->Branch("gv_nTkHi", gv_nTkHi, "gv_nTkHi[gv_n]/S");
  tree->Branch("gv_nTkLo", gv_nTkLo, "gv_nTkLo[gv_n]/S");
}

bool GlobeGenVertices::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GenParticleCollection> gpH;
  iEvent.getByLabel(genParticlesColl, gpH);

  gv_p3->Clear();
  gv_pos->Clear();
  gv_n = 0;

  const float lowPtThrGenVtx = 0.1;
  const float highPtThrGenVtx = 0.5;

  for(reco::GenParticleCollection::const_iterator it_gen = 
	gpH->begin(); it_gen!= gpH->end(); it_gen++){   
    if( it_gen->status() != 3 || !(it_gen->vx()!=0. || it_gen->vy()!=0. || it_gen->vx()!=0.)  ) { continue; }

    // check for duplicate vertex
    bool duplicate = false;
    for(Int_t itv = 0; itv < gv_n; itv++) {
      TVector3 * checkVtx = (TVector3 *) gv_pos->At(itv);
      if( (fabs(it_gen->vx()-checkVtx->X())<1e-5) &&  (fabs(it_gen->vy()-checkVtx->Y())<1e-5) && (fabs(it_gen->vz()-checkVtx->Z())<1e-5)) {
	duplicate = true;
	break;
      }
    }

    if (duplicate) continue;
    
    new((*gv_pos)[gv_n]) TVector3();
    ((TVector3 *) gv_pos->At(gv_n))->SetXYZ(it_gen->vx(), it_gen->vy(), it_gen->vz());
    
    TVector3 * this_gv_pos = (TVector3 *) gv_pos->At(gv_n);
    TVector3 p3(0,0,0);
    
    gv_sumPtLo[gv_n] = 0;
    gv_nTkLo[gv_n] = 0;
    gv_sumPtHi[gv_n] = 0;
    gv_nTkHi[gv_n] = 0;

    for(reco::GenParticleCollection::const_iterator part = gpH->begin(); part!= gpH->end(); part++){   
      if( part->status() == 1 && part->charge() != 0 && fabs(part->eta())<2.5 &&
	  ( fabs(part->vx()-this_gv_pos->X())<1.e-5 && fabs(part->vy()-this_gv_pos->Y())<1.e-5 && fabs(part->vz()-this_gv_pos->Z())<1.e-5 ) )  {
	
	TVector3 m(part->px(),part->py(),part->pz());
	p3 += m;
	if( m.Pt() > lowPtThrGenVtx ) {
	  gv_sumPtLo[gv_n] += m.Pt();
	  gv_nTkLo[gv_n] += 1;
	  if( m.Pt() > highPtThrGenVtx ) {
	    gv_sumPtHi[gv_n] += m.Pt();
	    gv_nTkHi[gv_n] += 1;
	  }
	}
      }
    }
    new((*gv_p3)[gv_n]) TVector3();
    ((TVector3 *) gv_p3->At(gv_n))->SetXYZ(p3.X(),p3.Y(),p3.Z());

    gv_n++;
  }

  return true;
}

