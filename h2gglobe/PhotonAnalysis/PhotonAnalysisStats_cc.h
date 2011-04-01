
void LoopAll::myStatPhotonAnalysisRed(Util *ut, int jentry){

  if (PADEBUG)
      cout << "myStatPhotonAnalysRed START" << endl;

  counters[0]++;
	
  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_calopos->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  std::vector<Elec> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;
  for (int i=0; i<pho_n && i < 4; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());
    //PreSelection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.1
       && pho_trksumpthollowconedr04[i] < 3*(1.5 + 0.001*pt)
       && pho_ecalsumetconedr04[i] < 2*(2.0 + 0.006*pt)
       && pho_hcalsumetconedr04[i] < 2*(2.0 + 0.0025*pt)
       &&((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5))) 
       ) {
         Elec candidate;
         candidate.p4 		= p4;
	 candidate.calopos	= calopos;
         candidate.pixSeed 	= pho_haspixseed[i];
         candidate.trkIso 	= pho_trksumpthollowconedr04[i];
         candidate.ecalIso 	= pho_ecalsumetconedr04[i];
         candidate.hcalIso 	= pho_hcalsumetconedr04[i];
         candidate.sieie 	= pho_sieie[i];
         candidate.hoe 		= pho_hoe[i];
         preselected_photons.push_back(candidate);
       }
  }
  
  
  

  float best_mass = 0.;
  float best_pt ;

  //Event Selection
  int n_preselected_pho = preselected_photons.size();

  if (n_preselected_pho > 1 ){

     bool dr_lead_match=false;
     bool dr_nlead_match=false;

     Elec leading, nleading;

     if (preselected_photons[0].p4->Pt() >= preselected_photons[1].p4->Pt()){

             leading  = preselected_photons[0];
             nleading = preselected_photons[1];
     } else {

             leading  = preselected_photons[1];
             nleading = preselected_photons[0];
     }


       if (leading.p4->Pt() > 40.){
           // -------------------------------------------------------
           TLorentzVector Higgs = (*(preselected_photons[0].p4))
                                 +(*(preselected_photons[1].p4));
           float mass = Higgs.M();
           float h_pt = Higgs.Pt();
             if (mass > 100. && mass < 150.){
             //Good event, passes preselection and acceptance cuts

             int pass_selection[2];
             int pass_isolation[2];
    
            //Now do selection on leading photon
             pass_selection[0] = 
                   leading.hoe < 0.02
                && leading.ecalIso < (2.0 + 0.006*leading.p4->Pt())
                && leading.hcalIso < (2.0 + 0.0025*leading.p4->Pt())
                && (((leading.sieie < 0.01)  && (leading.calopos->Eta() < 1.4442)) 
                || (( leading.sieie < 0.028)
                   && ((fabs(leading.calopos->Eta()) < 2.5) && (fabs(leading.calopos->Eta()) > 1.566))) );
             pass_isolation[0] =  leading.trkIso < (1.5 + 0.001*leading.p4->Pt());
	 
             //Selection on next to leading photon
             pass_selection[1] = 
                   (nleading.hoe < 0.02)
                && (nleading.ecalIso < (2.0 + 0.006*nleading.p4->Pt()))
                && (nleading.hcalIso < (2.0 + 0.0025*nleading.p4->Pt()))
                && (((nleading.sieie < 0.01)  && (nleading.calopos->Eta() < 1.4442)) 
                || (( nleading.sieie < 0.028)
                  && ((fabs(nleading.calopos->Eta()) < 2.5) && (fabs(nleading.calopos->Eta()) > 1.566))) );
             pass_isolation[1] =  nleading.trkIso < (1.5 + 0.001*nleading.p4->Pt());

           
             if (pass_selection[0] && pass_isolation[0]){
               if (pass_selection[1] && pass_isolation[1]){
                 best_mass = mass;
               }
             }


           }
     }
   }

     rooContainer->SetRealVar("mass",best_mass);
     rooContainer->SetRealVar("mass_sig",best_mass);
  
 
   if(PADEBUG)
     cout<<"myStatPhotonAnalysisRed END"<<endl;

}
