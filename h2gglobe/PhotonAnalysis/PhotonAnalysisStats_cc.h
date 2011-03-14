
void LoopAll::myStatPhotonAnalysisRed(Util *ut, int jentry){

  if (PADEBUG)
      cout << "myStatPhotonAnalysRed START" << endl;

  counters[0]++;
	
  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  std::vector<TLorentzVector *> selected_photons;  

  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    float pt  = p4->Pt();
    float eta = fabs(p4->Eta());
    if(
	(! (pho_haspixseed[i]==1))
     &&	(pt > 30. )
     && (pho_hoe[i] < 0.02)
     && (pho_trksumpthollowconedr04[i] < (1.5 + 0.001*pt))
     && (pho_ecalsumetconedr04[i] < (2.0 + 0.006*pt)	)
     && (pho_hcalsumetconedr04[i] < (2.0 + 0.0025*pt)	)
     && (   ((pho_sieie[i] < 0.01) && (eta < 1.44)) 
	 || ((pho_sieie[i] <0.028) && ((eta < 2.5) && (eta > 1.57))) 
	)
     && ((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5)) )
 
	){
		selected_photons.push_back(p4);
     }
  }

  int n_selected_pho = selected_photons.size();

  // now make sure there is one photon with Pt > 40 in the candidates
  float best_mass = 0.;
  
  if (n_selected_pho > 1) {
	
	for (int i=0; i< n_selected_pho-1; i++)
	  for (int j=i+1;j<n_selected_pho;j++){
		
		if (selected_photons[i]->Pt() > 40. 
		 || selected_photons[j]->Pt() > 40.){
				
		  TLorentzVector Higgs = (*selected_photons[i])+(*selected_photons[j]);
		  float mass = Higgs.M();
		  if (mass > best_mass) best_mass = mass;
		}
	  }	
  }
     rooContainer->SetRealVar("mass",best_mass);
  
 
   if(PADEBUG)
     cout<<"myStatPhotonAnalysisRed END"<<endl;

}
