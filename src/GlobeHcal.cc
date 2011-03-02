#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHcal.h"


GlobeHcal::GlobeHcal(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  hcalBEColl =  iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");
  hcalFColl =  iConfig.getParameter<edm::InputTag>("HcalHitsFColl");
  hcalHoColl =  iConfig.getParameter<edm::InputTag>("HcalHitsHoColl");
  
  doHFHcal = iConfig.getParameter<bool>("doHFHcal");

  debug_level = iConfig.getParameter<int>("Debug_Level");

  edm::ParameterSet psetHcal = iConfig.getParameter<edm::ParameterSet>("HcalHitsCuts");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeHcal::defineBranch(TTree* tree) {

  // think about changing branch names for duplicate collections
  hc_p4 = new TClonesArray("TLorentzVector", MAX_HCALHITS);
  tree->Branch("hc_p4", "TClonesArray", &hc_p4, 32000, 0);
  tree->Branch("hc_n", &hc_n, "hc_n/I");
  tree->Branch("hc_type", &hc_type, "hc_type[hc_n]/I");

}

bool GlobeHcal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,
			GlobeLeptons *lep, GlobeElectrons *el, GlobeMuons *mu,
			GlobePhotons *pho, GlobeTracks *tk) {

  edm::Handle<HBHERecHitCollection> hcalBEH;
  edm::Handle<HFRecHitCollection> hcalFH;
  edm::Handle<HORecHitCollection> hcalHoH;

  iEvent.getByLabel(hcalBEColl, hcalBEH);
  iEvent.getByLabel(hcalFColl, hcalFH);
  iEvent.getByLabel(hcalHoColl, hcalHoH);

  edm::ESHandle<CaloGeometry> pG;
  const CaloGeometry* caloGeom;

  iSetup.get<CaloGeometryRecord>().get(pG);

  caloGeom = pG.product();

  if (debug_level > 9) {
    std::cout<<"GlobeHcal: hcalBEH->size() "<< hcalBEH->size()<<std::endl;
    std::cout<<"GlobeHcal: hcalFH->size() "<< hcalFH->size()<<std::endl;
    std::cout<<"GlobeHcal: hcalHoH->size() "<< hcalHoH->size()<<std::endl;
  }
  // now have collections

  hc_p4->Clear();

  hc_n = 0;
  Int_t hcB_n = 0;
  Int_t hcF_n = 0;
  Int_t hcHo_n = 0;

  TVector3 vtxPos(0,0,0);
  TVector3 vtxMom(1,1,0);
  TVector3 calPos(0,0,0);

  for(unsigned int i=0; i<hcalBEH->size(); i++) {

    if (hcB_n >= MAX_HCALHITS) {
      std::cout << "GlobeHcal: 1 WARNING TOO MANY HCAL HITS: now " 
		<< hcB_n + hcF_n + hcHo_n + 1<< " (allowed " 
		<< MAX_HCALHITS << ")" << std::endl;
      break;
    }

    HBHERecHitRef hc(hcalBEH, i);
    if(i<10) {
      if(debug_level == -17) {
	std::cout << "Hcal HBHE "<<i<<" "<<hc->energy()<<std::endl;
      }
    }


    GlobalPoint posi = caloGeom->getPosition(hc->id());
    calPos.SetXYZ(posi.x(),posi.y(),posi.z());

    for(int j=0;j<lep->lpt_n;++j) { //Begin Lepton Loop List 
	    

      //The ordering should not be used ... it is confusing Marco
      //What does this mean??? Marco 
      //if(j == lep->lpt_emu_n) 
      //j=lep->lpt_el_n+lep->lpt_mu_n; //Skip bad electrons
      //What does this mean??? Marco 
      //if (j >= lep->lpt_n)
      //continue;

      //Why this? If we save a lepton we should also keep the hits...
      //if( ((TLorentzVector*)lep->lpt_p4->At(j))->Pt() < 5 ) continue; //Skip Low Pt Leptons

      if(abs(lep->lpt_pdgid[j]) == 11) {  //accessing the electron collection                  
	vtxPos = *((TVector3*)( el->el_posvtx->At(lep->lpt_ind[j]) ));
	vtxMom = *((TVector3*)( el->el_momvtx->At(lep->lpt_ind[j]) ));
      } else if(abs(lep->lpt_pdgid[j]) == 13) {
	vtxPos = *((TVector3*)( mu->mu_posvtx->At(lep->lpt_ind[j]) ));
	vtxMom = *((TVector3*)( mu->mu_momvtx->At(lep->lpt_ind[j]) ));
      } else if(lep->lpt_pdgid[j] == 22) {
	vtxPos.SetXYZ(0,0,0); 
	//this may also be wrong... If the photon has a vertex, it should take the photon vertex
	vtxMom = ((TLorentzVector*)( pho->pho_p4->At(lep->lpt_ind[j] )))->Vect();
      } else {
	std::cout << "Error, Lepton is not a photon, electron, or muon!!!" << std::endl;
	continue;
      }
      if(abs(lep->lpt_pdgid[j]) == 13) {
	if(debug_level == -17) {
	  //std::cout << "Vtx Pos: (" << vtxPos.X() << "," << vtxPos.Y() << "," << vtxPos.Z() << ")" << std::endl;
	  //std::cout << "Vtx Mom: (" << vtxMom.X() << "," << vtxMom.Y() << "," << vtxMom.Z() << ")" << std::endl;
	}
      }

      if( gCUT->cut(*hc,vtxMom.DeltaR(calPos-vtxPos)) ) continue;
      else { // passes all cuts
	if(debug_level >1000) std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;
    
	new ((*hc_p4)[hcB_n]) TLorentzVector();
	((TLorentzVector *)hc_p4->At(hcB_n))
	  ->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());
	//->SetPtEtaPhiE(hc->energy()*sin(posi.theta()),
	//               posi.eta(),
	//               posi.phi(),
	//               hc->energy());

	hc_type[hcB_n] = 1;
	hcB_n++;
	break; //break out of lepton loop, already passed 
		
      } //End Passes All Cuts

    } //End Lepton List loop 

    if(lep->lpt_n==0) {
      if( !gCUT->cut(*hc,100.) ) {
	new ((*hc_p4)[hcB_n]) TLorentzVector(); 
	((TLorentzVector *)hc_p4->At(hcB_n))
	  ->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());
	//->SetPtEtaPhiE(hc->energy()*sin(posi.theta()),
	//               posi.eta(),
	//               posi.phi(),
	//               hc->energy());
		
	hc_type[hcB_n] = 1;
	hcB_n++;
      }
    }

  }


  if (doHFHcal) {

    for(unsigned int i=0; i<hcalFH->size(); i++) {

      if (hcB_n + hcF_n >= MAX_HCALHITS) {
	std::cout << "GlobeHcal: 2 WARNING TOO MANY HCAL HITS: " << hcB_n + hcF_n + hcHo_n + 1 << " (allowed " << MAX_HCALHITS << ")" << std::endl;
	break;
      }

      HFRecHitRef hc(hcalFH, i);
      if(i<10) {
	if(debug_level == -17) {
	  std::cout << "Hcal HF   "<<i<<" "<<hc->energy()<<std::endl;
	}
      }
      // make the cuts
      if(gCUT->cut(*hc))continue; 
      // passed cuts

      if(debug_level >1000) std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;

      /*
	float phi =  caloGeom->getPosition(hc->id()).phi();
	float theta =  caloGeom->getPosition(hc->id()).theta();
	float en = hc->energy();
	float px = en*sin(theta)*cos(phi);
	float py = en*sin(theta)*sin(phi);
	float pz = en*cos(theta);
	new ((*hc_p4)[hcB_n + hcF_n]) TLorentzVector();
	((TLorentzVector *)hc_p4->At(hcB_n + hcF_n))->SetXYZT(px, py, pz, en);
      */

      GlobalPoint posi = caloGeom->getPosition(hc->id());
	    
      new ((*hc_p4)[hcB_n + hcF_n]) TLorentzVector();
      ((TLorentzVector *)hc_p4->At(hcB_n + hcF_n))->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());

      hc_type[hcB_n + hcF_n] = 2;
      hcF_n++;   
    }  
  }

  for(unsigned int i=0; i<hcalHoH->size(); i++) {
      
    if (hcB_n + hcF_n + hcHo_n >= MAX_HCALHITS) {
      std::cout << "GlobeHcal: WARNING TOO MANY HCAL HITS: now " << hcB_n + hcF_n + hcHo_n + 1<< " (allowed " << MAX_HCALHITS << ")" << std::endl;
      continue;
    }

    HORecHitRef hc(hcalHoH, i);
    if(i<10) {
      if(debug_level == -17) {
	std::cout << "Hcal HO   "<<i<<" "<<hc->energy()<<std::endl;
      }
    }
    // make the cuts
    if(gCUT->cut(*hc))continue; 
    // passed cuts

    if(debug_level >1000) std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;
    // passed cuts
    /*
      float phi =  caloGeom->getPosition(hc->id()).phi();
      float theta =  caloGeom->getPosition(hc->id()).theta();
      float en = hc->energy();
      float px = en*sin(theta)*cos(phi);
      float py = en*sin(theta)*sin(phi);
      float pz = en*cos(theta);

      new ((*hc_p4)[hcB_n + hcF_n + hcHo_n]) TLorentzVector();
      ((TLorentzVector *)hc_p4->At(hcB_n + hcF_n + hcHo_n))->SetXYZT(px, py, pz, en);
    */

    GlobalPoint posi = caloGeom->getPosition(hc->id());
        
    new ((*hc_p4)[hcB_n + hcF_n + hcHo_n]) TLorentzVector();
    ((TLorentzVector *)hc_p4->At(hcB_n + hcF_n + hcHo_n))->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());
    hc_type[hcB_n + hcF_n + hcHo_n] = 3;

    hcHo_n++;
  }

  hc_n = hcB_n + hcF_n + hcHo_n;

  //if (hc_n == 0)
  //  return false;

  return true;
}
