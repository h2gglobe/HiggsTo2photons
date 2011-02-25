#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePAT.h"




const double HIGH_ETA = 5.0;
const double MID_ETA = 3.5;
const double LOW_ETA = 2.5;

GlobePAT::GlobePAT(const edm::ParameterSet& iConfig):cuts_(iConfig) {
  
  eleLabel_ =iConfig.getUntrackedParameter<edm::InputTag>("electronTag");
  phoLabel_ =iConfig.getUntrackedParameter<edm::InputTag>("photonTag");
  jetLabel_ =iConfig.getUntrackedParameter<edm::InputTag>("jetTag");
  if(debug_level > 99 )
    std::cout << "GlobePAT: constructor" << std::endl;
  //  rechitbtag_ =iConfig.getUntrackedParameter<edm::InputTag>("rechitBTag");
  //rechitetag_ =iConfig.getUntrackedParameter<edm::InputTag>("rechitETag");
}

void GlobePAT::defineBranch(TTree* tree) {

  //set up TClonesArrays
  pat_el_p4 = new TClonesArray("TLorentzVector", MAX_PAT);
  pat_el_matchp4 = new TClonesArray("TLorentzVector", MAX_PAT);
  pat_el_motherp4 = new TClonesArray("TLorentzVector", MAX_PAT);
 
  pat_pho_p4 = new TClonesArray("TLorentzVector", MAX_PAT);
  pat_pho_matchp4 = new TClonesArray("TLorentzVector", MAX_PAT);
  pat_pho_motherp4 = new TClonesArray("TLorentzVector", MAX_PAT);

  pat_jet_p4 = new TClonesArray("TLorentzVector", MAX_PAT);
 
  
  //pat electron branches
  tree->Branch("pat_el_n", &pat_el_n, "pat_el_n/I");
  tree->Branch("pat_el_recoind", &pat_el_recoind, "pat_el_recoind[pat_el_n]/I");
  tree->Branch("pat_el_ismatched", &pat_el_ismatched, "pat_el_ismatched[pat_el_n]/I");
  tree->Branch("pat_el_isHEEP", &pat_el_isHEEP, "pat_el_isHEEP[pat_el_n]/I");
  tree->Branch("pat_el_HEEPbit", &pat_el_HEEPbit, "pat_el_HEEPbit[pat_el_n]/I");
  tree->Branch("pat_el_isHEEP_noDETA", &pat_el_isHEEP_noDETA, "pat_el_isHEEP_noDETA[pat_el_n]/I");
  tree->Branch("pat_el_p4","TClonesArray",&pat_el_p4,32000,0);
  tree->Branch("pat_el_matchp4","TClonesArray",&pat_el_matchp4,32000,0);
  tree->Branch("pat_el_motherp4","TClonesArray",&pat_el_motherp4,32000,0);
  tree->Branch("pat_el_motherpdgid", &pat_el_motherpdgid, "pat_el_motherpdgid[pat_el_n]/I");
  
  
  //pat photon branches
  tree->Branch("pat_pho_n", &pat_pho_n, "pat_pho_n/I");
  tree->Branch("pat_pho_recoind", &pat_pho_recoind, "pat_pho_recoind[pat_pho_n]/I");
  tree->Branch("pat_pho_ismatched", &pat_pho_ismatched, "pat_pho_ismatched[pat_pho_n]/I");
  tree->Branch("pat_pho_p4","TClonesArray",&pat_pho_p4,32000,0);
  tree->Branch("pat_pho_matchp4","TClonesArray",&pat_pho_matchp4,32000,0);
  tree->Branch("pat_pho_motherp4","TClonesArray",&pat_pho_motherp4,32000,0);
  tree->Branch("pat_pho_motherpdgid", &pat_pho_motherpdgid, "pat_pho_motherpdgid[pat_pho_n]/I");
  tree->Branch("pat_pho_TightID", &pat_pho_TightID, "pat_pho_TightID[pat_pho_n]/I");
  tree->Branch("pat_pho_TightID_EarlyData", &pat_pho_TightID_EarlyData, "pat_pho_TightID_EarlyData[pat_pho_n]/I");
  tree->Branch("pat_pho_LooseID", &pat_pho_LooseID, "pat_pho_LooseID[pat_pho_n]/I");
  tree->Branch("pat_pho_LooseID_EarlyData", &pat_pho_LooseID_EarlyData, "pat_pho_LooseID_EarlyData[pat_pho_n]/I");
  tree->Branch("pat_pho_IDBit_EarlyData", &pat_pho_IDBit_EarlyData, "pat_pho_IDBit_EarlyData[pat_pho_n]/I");
  

  //pat jet branches
  tree->Branch("pat_jet_n", &pat_jet_n, "pat_jet_n/I");
  tree->Branch("pat_jet_p4","TClonesArray",&pat_jet_p4,32000,0);
  tree->Branch("pat_jet_recoind", &pat_jet_recoind, "pat_jet_recoind[pat_jet_n]/I");
  if(debug_level > 99 )
    std::cout << "GlobePAT: defined branches" << std::endl;
}


bool GlobePAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeElectrons * ele, GlobePhotons * pho, GlobeJets * jet) {

  //Double_t dummy = -99.0;
  if(debug_level > 99 )
    std::cout << "GlobePAT: Start " << std::endl;

  pat_el_p4->Clear();
  pat_el_matchp4->Clear();
  pat_el_motherp4->Clear();
  pat_pho_p4->Clear();
  pat_pho_matchp4->Clear();
  pat_pho_motherp4->Clear();
  pat_jet_p4->Clear();
  pat_el_n = 0;
  pat_pho_n = 0;
  pat_jet_n = 0;
  
  //Do PAT electrons
  //Store some things, but also give index to el_std collection
  edm::Handle<edm::View<pat::Electron> > electronHandle;
  try{iEvent.getByLabel(eleLabel_,electronHandle);}
  catch(...){std::cout<<"Didn't find the pat::electron handle"<<std::endl;}
  const edm::View<pat::Electron> & electrons = *electronHandle;
  
  for(edm::View<pat::Electron>::const_iterator pat = electrons.begin(); pat!=electrons.end(); ++pat)
    {
      if(debug_level > 99 )
	std::cout << "Inside PAT electrons "<<"\t"<<pat_el_n << std::endl;
      new ((*pat_el_p4)[pat_el_n]) TLorentzVector();
      new ((*pat_el_matchp4)[pat_el_n]) TLorentzVector();
      new ((*pat_el_motherp4)[pat_el_n]) TLorentzVector();
     
      ((TLorentzVector *)pat_el_p4->At(pat_el_n))->SetXYZT(pat->px(), pat->py(), pat->pz(), pat->energy());
      
      TLorentzVector * elep4 = (TLorentzVector *) pat_el_p4->At(pat_el_n);
      pat_el_recoind[pat_el_n] = globeMatchPATtoReco(ele,elep4);
     
      //do some gen matching
      /*
      if(pat->genParticleRef().isNonnull()){
	((TLorentzVector *)pat_el_matchp4->At(pat_el_n))->SetXYZT(pat->genLepton()->px(),pat->genLepton()->py(),pat->genLepton()->pz(),pat->genLepton()->energy());
	pat_el_ismatched[pat_el_n] = 1;
      
	if(pat->genLepton()->mother() !=NULL){
	  reco::Candidate * mymother;
	  if(pat->genLepton()->mother()->status() != 3 )  mymother = (reco::Candidate *)pat->genLepton()->mother();
	  else  mymother = (reco::Candidate *)(pat->genLepton()->mother())->mother();  
	  
	  ((TLorentzVector *)pat_el_motherp4->At(pat_el_n))->SetXYZT(mymother->px(),mymother->py(),mymother->pz(),mymother->energy());
	
	  pat_el_motherpdgid[pat_el_n] = mymother->pdgId();
	}
	else{
	  pat_el_motherpdgid[pat_el_n] = -99;
	  ((TLorentzVector *)pat_el_motherp4->At(pat_el_n))->SetXYZT(-99.0,-99.0,-99.0,-99.0);
	}
	
      }
            
      else {
	((TLorentzVector *)pat_el_matchp4->At(pat_el_n))->SetXYZT(-99.0,-99.0,-99.0,-99.0);
	pat_el_ismatched[pat_el_n] = 0;
      }   
      */
      pat_el_ismatched[pat_el_n] = 0;
         //heep::CutCodes::passCuts(heepBitResult,~heep::CutCodes::DETAIN);
      //try to check if passes PAT
      //   std::cout<<"getcutcode  "<<cuts_.getCutCode(*pat)<<"\t" <<et<< std::endl;
      pat_el_HEEPbit[pat_el_n] = cuts_.getCutCode(*pat);
    
      if(fabs(pat->eta()) > 1.442 ){
        
	if( cuts_.passCuts(*pat,~heep::CutCodes::DETAIN) )
	  pat_el_isHEEP_noDETA[pat_el_n] =1 ;
	else 
	   pat_el_isHEEP_noDETA[pat_el_n] =0 ;
	
	if(cuts_.passCuts(*pat))
	  pat_el_isHEEP[pat_el_n] =1;
	else
	  pat_el_isHEEP[pat_el_n] =0;
	
      }
      else if (fabs(pat->eta()) < 1.442 &&  cuts_.passCuts(*pat)) {
	pat_el_isHEEP[pat_el_n] =1;
	pat_el_isHEEP_noDETA[pat_el_n] =1 ;
	
      }
      else if (fabs(pat->eta()) < 1.442) {
	pat_el_isHEEP[pat_el_n] =0;	
	pat_el_isHEEP_noDETA[pat_el_n] =0;	
      }
      /*   if( pat_el_isHEEP[pat_el_n] ==1 && cuts_.getCutCode(*pat) != 0 ){
       
       std::cout<<"Tight IDs are different = "<<pat_el_isHEEP[pat_el_n]<<"\t"<<cuts_.getCutCode(*pat)<<std::endl;

      }
      */
      pat_el_n++;
    }
  
  //    for(size_t eleNr=0;eleNr<eles.size();eleNr++){
  //  if(cuts_.passCuts(eles[eleNr])) nrPass_++;
  //  else nrFail_++;
  // }

   
   //Do PAT electrons
   //Store some things, but also give index to el_std collection
   edm::Handle<edm::View<pat::Photon> > phoHandle;
   try{iEvent.getByLabel(phoLabel_,phoHandle);}
   catch(...){std::cout<<"Didn't find the pat::photon handle"<<std::endl;}
   const edm::View<pat::Photon> & photons = *phoHandle;   

   for(edm::View<pat::Photon>::const_iterator pat = photons.begin(); pat!=photons.end(); ++pat) {
      if(debug_level > 99 )
	std::cout << "Inside PAT photons "<<"\t"<<pat_pho_n << std::endl;
    
     new ((*pat_pho_p4)[pat_pho_n]) TLorentzVector();
     new ((*pat_pho_matchp4)[pat_pho_n]) TLorentzVector();
     new ((*pat_pho_motherp4)[pat_pho_n]) TLorentzVector();
     ((TLorentzVector *)pat_pho_p4->At(pat_pho_n))->SetXYZT(pat->px(), pat->py(), pat->pz(), pat->energy());
     
     //       mypatPhoton_container.push_back(mypatPhoton);
     TLorentzVector * phop4 = (TLorentzVector *) pat_pho_p4->At(pat_pho_n);
     pat_pho_recoind[pat_pho_n] = globeMatchPATtoReco(pho,phop4);
     
     //Photon ID
     int phoIDBit = 0x0;

     if(pat->photonID("PhotonCutBasedIDTight") )  
       pat_pho_TightID[pat_pho_n] = 1;
     else   pat_pho_TightID[pat_pho_n] = 0;

     if(pat->photonID("PhotonCutBasedIDLoose") )  
       pat_pho_LooseID[pat_pho_n] = 1;
     else   pat_pho_LooseID[pat_pho_n] = 0;


     //early data
     if(pat->sigmaIetaIeta() > 0.013 ) phoIDBit |= 0x0001;
     if( pat->ecalRecHitSumEtConeDR04() > (4.2 + 0.004 * pat->et() ) ) phoIDBit |= 0x0002;
     if(pat->hcalTowerSumEtConeDR04() > (2.2 + 0.001 * pat->et() ) ) phoIDBit |= 0x0004;
     if(pat->trkSumPtHollowConeDR04() > (2.0 + 0.001 * pat->et() ) ) phoIDBit |= 0x0008;
     if( pat->hadronicOverEm() > 0.05 ) phoIDBit |= 0x0010;

     if(phoIDBit ==0 || phoIDBit ==1 ){
       pat_pho_LooseID_EarlyData[pat_pho_n] = 1;
     }
     else pat_pho_LooseID_EarlyData[pat_pho_n] = 0;

     if(phoIDBit ==0 ){
       pat_pho_TightID_EarlyData[pat_pho_n] = 1;
     }
     else pat_pho_TightID_EarlyData[pat_pho_n] = 0;

     
     pat_pho_IDBit_EarlyData[pat_pho_n] = phoIDBit;


     if( pat_pho_TightID_EarlyData[pat_pho_n] ==1 && phoIDBit !=0 ){
       
      std::cout<<"Tight IDs are different = "<<pat_pho_TightID_EarlyData[pat_pho_n]<<"\t"<<phoIDBit<<std::endl;

     }
     //do some gen matching
     /*
     if(pat->genParticleRef().isNonnull()){
       ((TLorentzVector *)pat_pho_matchp4->At(pat_pho_n))->SetXYZT(pat->genPhoton()->px(),pat->genPhoton()->py(),pat->genPhoton()->pz(),pat->genPhoton()->energy());
       pat_pho_ismatched[pat_pho_n] = 1;
       
       if(pat->genPhoton()->mother() !=NULL){
	 reco::Candidate * mymother;
	 if(pat->genPhoton()->mother()->status() != 3 )  mymother = (reco::Candidate *)pat->genPhoton()->mother();
	 else mymother = (reco::Candidate *)(pat->genPhoton()->mother())->mother();
	 
	 ((TLorentzVector *)pat_pho_motherp4->At(pat_pho_n))->SetXYZT(mymother->px(),mymother->py(),mymother->pz(),mymother->energy());
	 pat_pho_motherpdgid[pat_pho_n] = mymother->pdgId();
       }
       else{
	 pat_pho_motherpdgid[pat_pho_n] = -99;
	 ((TLorentzVector *)pat_pho_motherp4->At(pat_pho_n))->SetXYZT(-99.0,-99.0,-99.0,-99.0);
	 }
     }
     else {
       ((TLorentzVector *)pat_pho_matchp4->At(pat_pho_n))->SetXYZT(-99.0,-99.0,-99.0,-99.0);
       pat_pho_ismatched[pat_pho_n] = 0;
     }
     */
     pat_pho_ismatched[pat_pho_n] = 0;
     pat_pho_n++;
   }
   

   //do pat jets
   edm::Handle<edm::View<pat::Jet> > jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);
   const edm::View<pat::Jet> & jets = *jetHandle;
   
   for(edm::View<pat::Jet>::const_iterator pat = jets.begin(); pat!=jets.end(); ++pat){
     if(debug_level > 99 )
       std::cout << "Inside PAT jets "<<"\t"<<pat_jet_n << std::endl;
     new ((*pat_jet_p4)[pat_jet_n]) TLorentzVector();
     ((TLorentzVector *)pat_jet_p4->At(pat_jet_n))->SetXYZT(pat->px(), pat->py(), pat->pz(), pat->energy());
    
     TLorentzVector * jetp4 = (TLorentzVector *) pat_jet_p4->At(pat_jet_n); 
     pat_jet_recoind[pat_jet_n] = globeMatchPATtoReco(jet,jetp4); 
     

     
     pat_jet_n++;
   }

   /*
   Handle<EcalRecHitCollection> Brechit;//barrel
   Handle<EcalRecHitCollection> Erechit;//endcap
   
   try{ iEvent.getByLabel(rechitbtag_,Brechit); } catch(...) { std::cout<<"No handle to rechitEB found"<<std::endl;}
   try{ iEvent.getByLabel(rechitetag_,Erechit); } catch(...) { std::cout<<"No handle to rechitEE found"<<std::endl;}
   
   */
   


   //debug pring
   if (debug_level > 9)
     std::cout << "GlobeMet: PAT electron collection size: "<< electronHandle->size() << std::endl;
   
   if(debug_level > 99 )
     std::cout << "GlobePAT: End " << std::endl;
   
   return true;
   
}


int GlobePAT::globeMatchPATtoReco(GlobeElectrons * ele, TLorentzVector* p4)
{

  for (int i = 0; i< ele->el_n; i++){
    TLorentzVector * tempp4 = (TLorentzVector *) ele->el_p4->At(i);
    if(p4->DeltaR(*tempp4) == 0 ) {
      // std::cout<<"I found the reco ele to go with my PAT ele"<<std::endl;
      return i;
    }
  }

  std::cout<<"Didn't find an ele match, that's strange"<<std::endl;
  return -1;

}

int GlobePAT::globeMatchPATtoReco(GlobePhotons * pho, TLorentzVector* p4)
{

  for (int i = 0; i< pho->pho_n; i++){
    TLorentzVector * tempp4 = (TLorentzVector *) pho->pho_p4->At(i);
    if(p4->DeltaR(*tempp4) == 0 ) {
      // std::cout<<"I found the reco pho to go with my PAT ele"<<std::endl;
      return i;
    }
  }

  std::cout<<"Didn't find a pho match, that's strange"<<std::endl;
  return -1;

}

int GlobePAT::globeMatchPATtoReco(GlobeJets * jet, TLorentzVector* p4)
{
  //  std::cout<<"jet n = " << jet->jet_n<<std::endl;
  for (int i = 0; i< jet->jet_n; i++){
    TLorentzVector * tempp4 = (TLorentzVector *) jet->jet_p4->At(i);
    if(p4->DeltaR(*tempp4) <= 0.001 ) {
    //if(p4->Et() == tempp4->Et()){
      // std::cout<<"I found the reco pho to go with my PAT ele"<<std::endl;
      return i;
    }
  }

  return -1;

}

