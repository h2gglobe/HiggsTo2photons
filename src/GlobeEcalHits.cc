#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalHits.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
 
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
 
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"


GlobeEcalHits::GlobeEcalHits(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  doPreshowerHits = iConfig.getParameter<bool>("doPreshowerHits");

  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");
  if(doPreshowerHits) 
    ecalHitESColl = iConfig.getParameter<edm::InputTag>("EcalHitESColl");
  
  
  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeEcalHits::defineBranch(TTree* tree) {
  
  tree->Branch("ecalhit_n", &ecalhit_n, "ecalhit_n/I");
  tree->Branch("ecalhit_type", &ecalhit_type, "ecalhit_type[ecalhit_n]/S");
  tree->Branch("ecalhit_flag", &ecalhit_flag, "ecalhit_flag[ecalhit_n]/S");
  tree->Branch("ecalhit_time", &ecalhit_time, "ecalhit_time[ecalhit_n]/F");
  tree->Branch("ecalhit_detid", &ecalhit_detid, "ecalhit_detid[ecalhit_n]/I");
  tree->Branch("ecalhit_ieta", &ecalhit_ieta, "ecalhit_ieta[ecalhit_n]/S");
  tree->Branch("ecalhit_iphi", &ecalhit_iphi, "ecalhit_iphi[ecalhit_n]/S");
  tree->Branch("ecalhit_ix", &ecalhit_ix, "ecalhit_ix[ecalhit_n]/S");
  tree->Branch("ecalhit_iy", &ecalhit_iy, "ecalhit_iy[ecalhit_n]/S");
  tree->Branch("ecalhit_zside", &ecalhit_zside, "ecalhit_zside[ecalhit_n]/S");

  ecalhit_p4 = new TClonesArray("TLorentzVector", MAX_ECALRECHITS);
  tree->Branch("ecalhit_p4", "TClonesArray", &ecalhit_p4, 32000, 0);
}

//bool GlobeEcalHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,
//                            GlobeLeptons *lep, GlobeElectrons *el, GlobeMuons *mu,
//                            GlobePhotons *pho) {
bool GlobeEcalHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,
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
  
  // geometry initialization
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);

  const CaloSubdetectorGeometry* EB = geometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  const CaloSubdetectorGeometry* EE = geometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
  const CaloSubdetectorGeometry* ES = 0;
  if(doPreshowerHits) 
    ES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);

  edm::Handle<EBRecHitCollection> pEBRecHitH;
  edm::Handle<EERecHitCollection> pEERecHitH;
  edm::Handle<ESRecHitCollection> pESRecHitH; 

  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  if(doPreshowerHits)
    iEvent.getByLabel(ecalHitESColl, pESRecHitH);

  if (debug_level > 9)
    std::cout << "GlobeEcalHits: Ecal Hits collection size: "
              <<  pEBRecHitH->size() + pEERecHitH->size() << std::endl;

  ecalhit_p4->Clear();
  ecalhit_n = 0;

  //Assume we have a lepton list that gives us a vtxPos and vtxMom
  TVector3 vtxPos(0,0,0);
  TVector3 vtxMom(1,1,0);
  TVector3 calPos(0,0,0);

  for (size_t i = 0; i < pEBRecHitH->size(); ++i) { // Begin Barrel RecHit Loop
      
    if (ecalhit_n >= MAX_ECALRECHITS) {
      std::cout << "GlobeEcalHits: WARNING TOO MANY RECHITS: " 
                << pEBRecHitH->size() + pEERecHitH->size() 
                << " (allowed " << MAX_ECALRECHITS << ")" << std::endl;
      break;
    }
      
    EcalRecHitRef rh(pEBRecHitH, i);
    DetId id = rh->detid(); 

    if(i<10) {
      if(debug_level == -17) {
        std::cout << "Ecal Barr "<<i<<" "<<rh->energy()<<std::endl;
      }
    }
        
    const CaloCellGeometry* this_cell = EB->getGeometry(id);
    if (this_cell) { //Being RecHit exists in geometry
      GlobalPoint posi = this_cell->getPosition();
      calPos.SetXYZ(posi.x(),posi.y(),posi.z());

      for (unsigned int j=0; j<leptons; j++) {
	vtxPos = *((TVector3*)lptPos->At(j));
	vtxMom = *((TVector3*)lptMom->At(j));
	
        if(gCUT->cut(*rh,0,vtxMom.DeltaR(calPos-vtxPos))) 
          continue;                   
        else {                                             //Passes all cuts
	  
          new ((*ecalhit_p4)[ecalhit_n]) TLorentzVector(); 
          ((TLorentzVector *)ecalhit_p4->At(ecalhit_n))->SetXYZT(posi.x(),posi.y(),posi.z(),rh->energy());
	  
          ecalhit_type[ecalhit_n] = 0; 
          ecalhit_time[ecalhit_n] = rh->time();
          ecalhit_flag[ecalhit_n] = rh->recoFlag();
	  ecalhit_detid[ecalhit_n] = rh->detid();
	  ecalhit_ieta[ecalhit_n] = ((EBDetId)rh->detid()).ieta();
	  ecalhit_iphi[ecalhit_n] = ((EBDetId)rh->detid()).iphi();
	  ecalhit_zside[ecalhit_n] = ((EBDetId)rh->detid()).zside();
	  ecalhit_ix[ecalhit_n] = -9999;
	  ecalhit_iy[ecalhit_n] = -9999;	  
          ecalhit_n++;
          break; //break out of lepton loop, already passed 
        } //End Passes All Cuts
      } //End Lepton List loop

      if(leptons == 0) {
        if(!gCUT->cut(*rh,0,100.)) {
          new ((*ecalhit_p4)[ecalhit_n]) TLorentzVector(); 
          ((TLorentzVector *)ecalhit_p4->At(ecalhit_n))->SetXYZT(posi.x(),posi.y(),posi.z(),rh->energy());

          ecalhit_type[ecalhit_n] = 0;
          ecalhit_n++;  
	}
      }

    } //End RecHit Exists in geometry
  } // End barrel RecHit Loop

  for (size_t i = 0; i < pEERecHitH->size(); ++i) {
    if (ecalhit_n >= MAX_ECALRECHITS) {
      std::cout << "GlobeEcalHits: WARNING TOO MANY RECHITS: " 
                << pEBRecHitH->size() + pEERecHitH->size() 
                << " (allowed " << MAX_ECALRECHITS << ")" << std::endl;
      break;
    }
      
    EcalRecHitRef rh(pEERecHitH, i);
      
    DetId id = rh->detid();
    const CaloCellGeometry* this_cell = EE->getGeometry(id);
    if (this_cell) { //Being RecHit exists in geometry
      GlobalPoint posi = this_cell->getPosition();
      calPos.SetXYZ(posi.x(),posi.y(),posi.z());

      for (unsigned int j=0; j<leptons; j++) {
	vtxPos = *((TVector3*)lptPos->At(j));
	vtxMom = *((TVector3*)lptMom->At(j));
          
        if( gCUT->cut(*rh,1,vtxMom.DeltaR(calPos-vtxPos)) ) continue;
        else {                                             //Passes all cuts
          new ((*ecalhit_p4)[ecalhit_n]) TLorentzVector(); 
          ((TLorentzVector *)ecalhit_p4->At(ecalhit_n))->SetXYZT(posi.x(),posi.y(),posi.z(),rh->energy());
            
          ecalhit_type[ecalhit_n] = 1;
          ecalhit_time[ecalhit_n] = rh->time();
          ecalhit_flag[ecalhit_n] = (int)rh->recoFlag();
	  ecalhit_detid[ecalhit_n] = rh->detid();
	  ecalhit_ix[ecalhit_n] = ((EEDetId)rh->detid()).ix();
	  ecalhit_iy[ecalhit_n] = ((EEDetId)rh->detid()).iy();
	  ecalhit_zside[ecalhit_n] = ((EEDetId)rh->detid()).zside();
	  ecalhit_ieta[ecalhit_n] = -9999;
	  ecalhit_iphi[ecalhit_n] = -9999;
          ecalhit_n++;
          break; //break out of lepton loop, already passed 
            
        } //End Passes All Cuts
          
      } //End Lepton List loop 

      if(leptons == 0) {
        if(!gCUT->cut(*rh, 1, 100.)) {
          new ((*ecalhit_p4)[ecalhit_n]) TLorentzVector(); 
          ((TLorentzVector *)ecalhit_p4->At(ecalhit_n))->SetXYZT(posi.x(),posi.y(),posi.z(),rh->energy());

          ecalhit_type[ecalhit_n] = 2;
          ecalhit_n++;  
	}
      }
    }
  }
    
  if (doPreshowerHits) {
    for (size_t i = 0; i < pESRecHitH->size(); ++i) {
        
      if (ecalhit_n >= MAX_ECALRECHITS) {
        std::cout << "GlobeEcalHits: WARNING TOO MANY RECHITS: " 
                  << pEBRecHitH->size() + pEERecHitH->size()+pESRecHitH->size() 
                  << " (allowed " << MAX_ECALRECHITS << ")" << std::endl;
        break;
      }
        
      EcalRecHitRef rh(pESRecHitH, i);
        
      if(i<10) {
	if(debug_level == -17) {
	  std::cout << "Ecal Pres "<<i<<" "<<rh->energy()<<std::endl;
	}
      }

      DetId id = rh->detid();
      const CaloCellGeometry* this_cell = ES->getGeometry(id);
      if (this_cell) { //Being RecHit exists in geometry
        GlobalPoint posi = this_cell->getPosition();
          
        if( gCUT->cut(*rh,2,0) )
          continue;
        else {  //passes all cuts
          new ((*ecalhit_p4)[ecalhit_n]) TLorentzVector(); 
          ((TLorentzVector *)ecalhit_p4->At(ecalhit_n))->SetXYZT(posi.x(),posi.y(),posi.z(),rh->energy());
          //->SetPtEtaPhiE(rh->energy()*sin(posi.theta()),posi.eta(),posi.phi(),rh->energy());
            
          ecalhit_type[ecalhit_n] = 2; 
          ecalhit_time[ecalhit_n] = rh->time();
          ecalhit_flag[ecalhit_n] = rh->recoFlag();
	  ecalhit_detid[ecalhit_n] = rh->detid();
	  ecalhit_ieta[ecalhit_n] = -9999;
	  ecalhit_iphi[ecalhit_n] = -9999;
	  ecalhit_ix[ecalhit_n] = -9999;
	  ecalhit_iy[ecalhit_n] = -9999;
	  ecalhit_zside[ecalhit_n] = -9999;
          ecalhit_n++;
            
        } //End Passes All Cuts
      }
    }
  }

  delete lptPos;
  delete lptMom;

  return true;
}
