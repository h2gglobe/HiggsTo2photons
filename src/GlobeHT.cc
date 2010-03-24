#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHT.h"

const double HIGH_ETA = 5.0;
const double MID_ETA = 3.5;
const double LOW_ETA = 2.5;

GlobeHT::GlobeHT(const edm::ParameterSet& iConfig) {

  trackColl =  iConfig.getParameter<edm::InputTag>("TrackColl");

}

void GlobeHT::defineBranch(TTree* tree) {

  ht_trkvec = new TVector3(0.0,0.0,0.0);

  tree->Branch("ht_25", &ht_25, "ht_25/F");
  tree->Branch("ht_35", &ht_35, "ht_35/F");
  tree->Branch("ht_50", &ht_50, "ht_50/F");
  tree->Branch("ht_nomet25", &ht_nomet25, "ht_nomet25/F");
  tree->Branch("ht_nomet35", &ht_nomet35, "ht_nomet35/F");
  tree->Branch("ht_nomet50", &ht_nomet50, "ht_nomet50/F");
  tree->Branch("ht_trk", &ht_trk, "ht_trk/F");
  tree->Branch("ht_trkvec", "TVector3", &ht_trkvec, 32000, 0);

  tree->Branch("ht_2lpt_n", &ht_2lpt_n, "ht_2lpt_n/I");
  tree->Branch("ht_2lpt_inds", &ht_2lpt_inds, "ht_2lpt_inds[ht_2lpt_n][2]/I");
  tree->Branch("ht_4lpt_n", &ht_4lpt_n, "ht_4lpt_n/I");
  tree->Branch("ht_4lpt_inds", &ht_4lpt_inds, "ht_4lpt_inds[ht_4lpt_n][4]/I");

  tree->Branch("ht_2lpt25", &ht_2lpt25, "ht_2lpt25[ht_2lpt_n]/F");
  tree->Branch("ht_2lpt35", &ht_2lpt35, "ht_2lpt35[ht_2lpt_n]/F"); 
  tree->Branch("ht_2lpt50", &ht_2lpt50, "ht_2lpt50[ht_2lpt_n]/F");
  tree->Branch("ht_4lpt25", &ht_4lpt25, "ht_4lpt25[ht_4lpt_n]/F");
  tree->Branch("ht_4lpt35", &ht_4lpt35, "ht_4lpt35[ht_4lpt_n]/F");
  tree->Branch("ht_4lpt50", &ht_4lpt50, "ht_4lpt50[ht_4lpt_n]/F");
}

void GlobeHT::fillTrackHT(const edm::Event& iEvent) {

    edm::Handle<reco::TrackCollection> tkH; 
    iEvent.getByLabel(trackColl, tkH);

    ht_trk=0.;
    ht_trkvec->SetXYZ(0,0,0);

    for(reco::TrackCollection::size_type j = 0; j<tkH->size(); ++j) {
        reco::TrackRef track(tkH, j);

        TVector3 trkTemp(track->px(),track->py(),track->pz());
        ht_trk += track->pt();
        *ht_trkvec += trkTemp;
    }

}

void GlobeHT::fillCaloTowerHT(GlobeMET * theMET, GlobeCaloTowers * theCaloTowers) {

  ht_nomet25=0.;
  ht_nomet35=0.;
  ht_nomet50=0.;
  for(Int_t ict=0;ict!=theCaloTowers->ct_n;++ict) {
      TLorentzVector * ctp4= (TLorentzVector *) theCaloTowers->ct_p4->At(ict);
      float ctEta = fabs(ctp4->Eta());
      if(ctEta < LOW_ETA)ht_nomet25+=ctp4->Et();
      if(ctEta < MID_ETA)ht_nomet35+=ctp4->Et();
      if(ctEta < HIGH_ETA)ht_nomet50+=ctp4->Et();
  }
  ht_25=ht_nomet25 + theMET->met_met_mip;
  ht_35=ht_nomet35 + theMET->met_met_mip;
  ht_50=ht_nomet50 + theMET->met_met_mip;

}

void GlobeHT::fillLeptonHT(GlobeJets * theJets, GlobeMET * theMET, GlobeLeptons * theLeptons) {

  ht_2lpt_n = 0;
  ht_4lpt_n=0;

  if(!theJets)return;
  if(!theMET)return;
  if(!theLeptons)return;

  float minJetEt = 15.;
  float drjetsep = 0.5;

  if(theLeptons->lpt_emu_n > 1) {
    for( Int_t i=0;i!=theLeptons->lpt_emu_n; ++i) {
      TLorentzVector * ip4= (TLorentzVector *) theLeptons->lpt_p4->At(i);
      for( Int_t ii=i+1;ii!=theLeptons->lpt_emu_n; ++ii) {
        TLorentzVector * iip4= (TLorentzVector *) theLeptons->lpt_p4->At(ii);
        if(ht_2lpt_n == MAX_HT2-1) {
          std::cout << "MAX_HT2 too small!!  " << MAX_HT2 << std::endl;
          continue;
        }

        // have two leptons (i and ii)
	ht_2lpt25[ht_2lpt_n] = theMET->met_met_mip + ip4->Et() + iip4->Et();
	ht_2lpt35[ht_2lpt_n] = ht_2lpt25[ht_2lpt_n]; 	 
	ht_2lpt50[ht_2lpt_n] = ht_2lpt25[ht_2lpt_n];
        ht_2lpt_inds[ht_2lpt_n][0]=i;
        ht_2lpt_inds[ht_2lpt_n][1]=ii;

        for (int ijet = 0; ijet != theJets->jet_n; ++ijet){
          TLorentzVector * jetp4= (TLorentzVector *) theJets->jet_p4->At(ijet);
          // minimal requirements on jets
          if(fabs( jetp4->Eta() ) > HIGH_ETA )continue;
          if(jetp4->Et() > minJetEt) {
              // continue if jet is close to any leptons
              if( ip4->DeltaR(*jetp4) < drjetsep)continue;
              if( iip4->DeltaR(*jetp4) < drjetsep)continue;
          }
          // this jet is separated from the 2 leptons 
          ht_2lpt50[ht_2lpt_n]+=jetp4->Et();
          if(fabs(jetp4->Eta()) > MID_ETA )continue;
          ht_2lpt35[ht_2lpt_n]+=jetp4->Et();
          if(fabs(jetp4->Eta()) > LOW_ETA )continue;
          ht_2lpt25[ht_2lpt_n]+=jetp4->Et();
        }// end for-loop(ijet)
        //std::cout << "ht_2lpt_n: " << ht_2lpt_n << std::endl;
        ht_2lpt_n++;
      }// end for-loop(ii)
    }// end for-loop(i)


    if(theLeptons->lpt_emu_n > 3) {
      for( Int_t i=0;i!=theLeptons->lpt_emu_n; ++i) {
        TLorentzVector * ip4= (TLorentzVector *) theLeptons->lpt_p4->At(i);

        for( Int_t ii=i+1;ii!=theLeptons->lpt_emu_n; ++ii) {
          TLorentzVector * iip4= (TLorentzVector *) theLeptons->lpt_p4->At(ii);

          for( Int_t iii=ii+1;iii!=theLeptons->lpt_emu_n; ++iii) {
            TLorentzVector * iiip4= (TLorentzVector *) theLeptons->lpt_p4->At(iii);

            for( Int_t iiii=iii+1;iiii!=theLeptons->lpt_emu_n; ++iiii) {
              TLorentzVector * iiiip4= (TLorentzVector *) theLeptons->lpt_p4->At(iiii);

              if(ht_4lpt_n == MAX_HT4-1) {
                std::cout << "MAX_HT4 too small!!  " << MAX_HT4 << std::endl;
                continue;
              }

              // have four leptons (i ii iii and iiii)
              ht_4lpt25[ht_4lpt_n] = theMET->met_met_mip + ip4->Et() + iip4->Et() + iiip4->Et() + iiiip4->Et();
	      ht_4lpt35[ht_4lpt_n] = ht_4lpt25[ht_4lpt_n]; 	 
	      ht_4lpt50[ht_4lpt_n] = ht_4lpt25[ht_4lpt_n];
	      ht_4lpt_inds[ht_4lpt_n][0]=i;
              ht_4lpt_inds[ht_4lpt_n][1]=ii;
              ht_4lpt_inds[ht_4lpt_n][2]=iii;
              ht_4lpt_inds[ht_4lpt_n][3]=iiii;

              for (int ijet = 0; ijet != theJets->jet_n; ++ijet){
                TLorentzVector * jetp4= (TLorentzVector *) theJets->jet_p4->At(ijet);
                if(fabs( jetp4->Eta() ) > HIGH_ETA )continue;  // only consider jets below eta of 5.0
                if(jetp4->Et() > minJetEt) {
                    // continue if jet is close to any leptons
                    if( ip4->DeltaR(*jetp4) < drjetsep)continue;
                    if( iip4->DeltaR(*jetp4) < drjetsep)continue;
                    if( iiip4->DeltaR(*jetp4) < drjetsep)continue;
                    if( iiiip4->DeltaR(*jetp4) < drjetsep)continue;
                }
                // this jet is separated from the 4 leptons 
                ht_4lpt50[ht_4lpt_n]+=jetp4->Et();
                if(fabs( jetp4->Eta() ) > MID_ETA )continue;
                ht_4lpt35[ht_4lpt_n]+=jetp4->Et();
                if(fabs( jetp4->Eta() ) > LOW_ETA )continue;
                ht_4lpt25[ht_4lpt_n]+=jetp4->Et();
              }// end for-loop(ijet)
              //std::cout << "ht_4lpt_n: " << ht_4lpt_n << std::endl;
              ht_4lpt_n++;
            }// end for-loop(iiii)
          }// end for-loop(iii)
        }// end for-loop(ii)
      }// end for-loop(i)
    }// end if(lpt_emu_n > 3)
  }// end if(lpt_emu_n > 1)

}

