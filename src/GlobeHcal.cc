#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHcal.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

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

void GlobeHcal::defineBranch(GlobeAnalyzer* ana) {

  // think about changing branch names for duplicate collections
  hc_p4 = new TClonesArray("TLorentzVector", MAX_HCALHITS);
  ana->Branch("hc_p4", "TClonesArray", &hc_p4, 32000, 0);
  ana->Branch("hc_n", &hc_n, "hc_n/I");
  ana->Branch("hc_type", &hc_type, "hc_type[hc_n]/I");
}

bool GlobeHcal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,
			GlobeElectrons *el, GlobeMuons *mu, GlobePhotons *pho) {

  TClonesArray* lptPos = new TClonesArray("TVector3");
  TClonesArray* lptMom = new TClonesArray("TVector3");

  unsigned int leptons = 0;
  if (el) {
    for (int i=0; i<el->el_n; i++) {
      new((*lptPos)[leptons]) TVector3();
      ((TVector3*)lptPos->At(leptons))->SetXYZ(((TVector3*)el->el_posvtx->At(i))->X(), ((TVector3*)el->el_posvtx->At(i))->Y(), ((TVector3*)el->el_posvtx->At(i))->Z());
      new((*lptMom)[leptons]) TVector3();
      ((TVector3*)lptMom->At(leptons))->SetXYZ(((TVector3*)el->el_momvtx->At(i))->X(), ((TVector3*)el->el_momvtx->At(i))->Y(), ((TVector3*)el->el_momvtx->At(i))->Z());
      leptons++;
    }
  }

  if (mu) {
    for (int i=0; i<mu->mu_n; i++) {
      new((*lptPos)[leptons]) TVector3();
      ((TVector3*)lptPos->At(leptons))->SetXYZ(((TVector3*)mu->mu_posvtx->At(i))->X(), ((TVector3*)mu->mu_posvtx->At(i))->Y(), ((TVector3*)mu->mu_posvtx->At(i))->Z());
      new((*lptMom)[leptons]) TVector3();
      ((TVector3*)lptMom->At(leptons))->SetXYZ(((TVector3*)mu->mu_momvtx->At(i))->X(), ((TVector3*)mu->mu_momvtx->At(i))->Y(), ((TVector3*)mu->mu_momvtx->At(i))->Z());
      leptons++;
    }
  }

  if (pho) {
    for (int i=0; i<pho->pho_n; i++) {
      new((*lptPos)[leptons]) TVector3();
      ((TVector3*)lptPos->At(leptons))->SetXYZ(0, 0, 0);
      new((*lptMom)[leptons]) TVector3();
      TVector3  temp = (((TLorentzVector*)(pho->pho_p4->At(i)))->Vect());
      ((TVector3*)lptMom->At(leptons))->SetXYZ(temp.X(), temp.Y(), temp.Z());
      leptons++;
    }
  }

  
  edm::Handle<HBHERecHitCollection> hcalBEH;
  edm::Handle<HFRecHitCollection> hcalFH;
  edm::Handle<HORecHitCollection> hcalHoH;
  iEvent.getByLabel(hcalBEColl, hcalBEH);
  
  edm::ESHandle<CaloGeometry> pG;
  const CaloGeometry* caloGeom;
  
  iSetup.get<CaloGeometryRecord>().get(pG);
  
  caloGeom = pG.product();

  if (debug_level > 9) {
    std::cout<<"GlobeHcal: hcalBEH->size() "<< hcalBEH->size()<<std::endl;
    //std::cout<<"GlobeHcal: hcalFH->size() "<< hcalFH->size()<<std::endl;
    //std::cout<<"GlobeHcal: hcalHoH->size() "<< hcalHoH->size()<<std::endl;
  }

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

    GlobalPoint posi = caloGeom->getPosition(hc->id());
    calPos.SetXYZ(posi.x(),posi.y(),posi.z());

    for (unsigned int j=0; j<leptons; j++) {
      vtxPos = *((TVector3*)lptPos->At(j));
      vtxMom = *((TVector3*)lptMom->At(j));
      
      if( gCUT->cut(*hc,vtxMom.DeltaR(calPos-vtxPos))) 
	continue;
      else { // passes all cuts
	if(debug_level >1000) std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;
    
	new ((*hc_p4)[hcB_n]) TLorentzVector();
	((TLorentzVector *)hc_p4->At(hcB_n))->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());

	hc_type[hcB_n] = 1;
	hcB_n++;
	break; //break out of lepton loop, already passed 		
      } //End Passes All Cuts
    } //End Lepton List loop 

    if(leptons == 0) {
      if(!gCUT->cut(*hc,100.)) {
	new ((*hc_p4)[hcB_n]) TLorentzVector(); 
	((TLorentzVector *)hc_p4->At(hcB_n))->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());
		
	hc_type[hcB_n] = 1;
	hcB_n++;
      }
    }
  }

  if (doHFHcal) {
    iEvent.getByLabel(hcalFColl, hcalFH);
    iEvent.getByLabel(hcalHoColl, hcalHoH);

    for(unsigned int i=0; i<hcalFH->size(); i++) {

      if (hcB_n + hcF_n >= MAX_HCALHITS) {
	std::cout << "GlobeHcal: 2 WARNING TOO MANY HCAL HITS: " << hcB_n + hcF_n + hcHo_n + 1 << " (allowed " << MAX_HCALHITS << ")" << std::endl;
	break;
      }

      HFRecHitRef hc(hcalFH, i);

      if(gCUT->cut(*hc))
	continue; 
      if(debug_level >1000) 
	std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;

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
    if(gCUT->cut(*hc))
      continue; 
    if(debug_level >1000) 
      std::cout<<"GlobeHcal: energy "<<hc->energy()<<std::endl;

    GlobalPoint posi = caloGeom->getPosition(hc->id());
        
    new ((*hc_p4)[hcB_n + hcF_n + hcHo_n]) TLorentzVector();
    ((TLorentzVector *)hc_p4->At(hcB_n + hcF_n + hcHo_n))->SetXYZT(posi.x(),posi.y(),posi.z(),hc->energy());
    hc_type[hcB_n + hcF_n + hcHo_n] = 3;

    hcHo_n++;
  }

  hc_n = hcB_n + hcF_n + hcHo_n;


  delete lptPos;
  delete lptMom;

  return true;
}
