#define PADEBUG 0

void LoopAll::InitRealPhotonAnalysis(Util *ut, int typerun) {

  // Book histos only if not reduce step
  if (typerun != 1) {    
   for (int ind=0;ind<ut->ntypes;ind++){
    // Standard Kinematic Plots -- These have Full Selection (ie lead/sub-lead photon)
    // Have a total, and one for each category (currently 4 categories)
    histoContainer[ind]->Add("h_mass",";M #gamma #gamma - GeV/c^{2};", 10, 100, 150);
    histoContainer[ind]->Add("h_pt",";p_{T} #gamma #gamma - GeV/c;", 50, 0, 50);
    //histoContainer[ind]->Add("h_n_sel",";N pre-selected photons;", 6, 0, 6);

    //Category 1 - R9_min < 0.93, |eta_max| in EB
    histoContainer[ind]->Add("h_mass_cat1",";M #gamma #gamma - GeV/c^{2};", 10, 100, 150);
    histoContainer[ind]->Add("h_pt_cat1",";p_{T} #gamma #gamma - GeV/c;", 50, 0, 50);

    //Category 2 - R9_min > 0.93, |eta_max| in EB
    histoContainer[ind]->Add("h_mass_cat2",";M #gamma #gamma - GeV/c^{2};", 10, 100, 150);
    histoContainer[ind]->Add("h_pt_cat2",";p_{T} #gamma #gamma - GeV/c;", 50, 0, 50);

    //Category 3 - R9_min < 0.93, |eta_max| in EE
    histoContainer[ind]->Add("h_mass_cat3",";M #gamma #gamma - GeV/c^{2};", 10, 100, 150);
    histoContainer[ind]->Add("h_pt_cat3",";p_{T} #gamma #gamma - GeV/c;", 50, 0, 50);

    //Category 4 - R9_min > 0.93, |eta_max| in EE
    histoContainer[ind]->Add("h_mass_cat4",";M #gamma #gamma - GeV/c^{2};", 10, 100, 150);
    histoContainer[ind]->Add("h_pt_cat4",";p_{T} #gamma #gamma - GeV/c;", 50, 0, 50);
    // -------------------------------------------------------------------------------

    // These have simply 30 GeV pt and fiducial eta Cuts
    histoContainer[ind]->Add("h_pho_pt",";p_{T} ^{#gamma};", 70, 30, 100);
    histoContainer[ind]->Add("h_pho_eta",";|#eta| ^{#gamma};", 25, 0, 2.5);

    // Isolation / ID variables for N-1 Plots, (range is to 3 * cut value)
    // These will include all photons which pass 30 GeV and |eta| restriction
    histoContainer[ind]->Add("h_pho_hoe",";H/E;", 60, 0, 0.06);
    histoContainer[ind]->Add("h_pho_sieie_eb",";#sigma_{i #eta i #eta} EB;", 30, 0, 0.03);
    histoContainer[ind]->Add("h_pho_sieie_ee",";#sigma_{i #eta i #eta} EE;", 90, 0, 0.09);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04",";trk iso;", 45, 0, 4.5);
    histoContainer[ind]->Add("h_pho_ecalsumetconedr04",";ecal iso;", 60, 0, 6.0);
    histoContainer[ind]->Add("h_pho_hcalsumetconedr04",";hcal iso;", 60, 0, 6.0);

    // Isolation Histograms for Pass/Fail distributions ---------------------------------------------------
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_leading_pass","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_leading_fail","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_nleading_pass","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_nleading_fail","TrkIso hollow DR04", 20, 0, 3);
    // ----------------------------------------------------------------------------------------------------

    // Gen-Level plots
    histoContainer[ind]->Add("h_dr_lead",150,0,0.05);
    histoContainer[ind]->Add("h_dr_nlead",150,0,0.05);

    // Correlation Plots -----------------------------------------------------------//

    histoContainer[ind]->Add("h_corr_trckiso_vs_hoe",";trkiso;H/E",  50, 0., 3., 50,0.,0.04 );
    histoContainer[ind]->Add("h_corr_trckiso_vs_sieie",";trkiso;sieie",  50, 0., 3., 50,0.,0.06 );
    histoContainer[ind]->Add("h_corr_ecaliso_vs_hoe",";ecaliso;H/E",  50, 0., 3., 50,0.,0.04 );
    histoContainer[ind]->Add("h_corr_ecaliso_vs_sieie",";ecaliso;sieie",  50, 0., 3., 50,0.,0.06 );
    histoContainer[ind]->Add("h_corr_hcaliso_vs_hoe",";hcaliso;H/E",  50, 0., 3., 50,0.,0.04 );
    histoContainer[ind]->Add("h_corr_hcaliso_vs_sieie",";hcaliso;sieie",  50, 0., 3., 50,0.,0.06 );
    
 
    // double-sideband method histograms -> All
    histoContainer[ind]->Add("h_selected_asymm", 2,0,2);
    histoContainer[ind]->Add("h_sideband_leading",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading", 2, -0.5, 1.5, 2, -0.5, 1.5);
 
    //Category 1 - R9_min < 0.93, |eta_max| in EB
    histoContainer[ind]->Add("h_selected_asymm_cat1", 2,0,2);
    histoContainer[ind]->Add("h_sideband_leading_cat1",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading_cat1", 2, -0.5, 1.5, 2, -0.5, 1.5);

    //Category 2 - R9_min > 0.93, |eta_max| in EB
    histoContainer[ind]->Add("h_selected_asymm_cat2", 2,0,2);
    histoContainer[ind]->Add("h_sideband_leading_cat2",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading_cat2", 2, -0.5, 1.5, 2, -0.5, 1.5);

    //Category 3 - R9_min > 0.93, |eta_max| in EE
    histoContainer[ind]->Add("h_selected_asymm_cat3", 2,0,2);
    histoContainer[ind]->Add("h_sideband_leading_cat3",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading_cat3", 2, -0.5, 1.5, 2, -0.5, 1.5);

    //Category 4 - R9_min < 0.93, |eta_max| in EE
    histoContainer[ind]->Add("h_selected_asymm_cat4", 2,0,2);
    histoContainer[ind]->Add("h_sideband_leading_cat4",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading_cat4", 2, -0.5, 1.5, 2, -0.5, 1.5);

   }
  }

    if (typerun == 3){
 
       rooContainer->AddRealVar("mass",100.,200.);
       rooContainer->AddRealVar("mass_sig",115.,125.);
       rooContainer->AddRealVar("mu1",-0.04,-1.0,-0.001);
       rooContainer->AddRealVar("mu2",-0.04,-1.0,-0.001);
       rooContainer->AddRealVar("mu3",-0.04,-1.0,-0.001);
 
       // make the vector of parameters as dependants of
       // the function to be fitted
       // ----------------------------------------------//
       std::vector<const char*> pars(2,"t");
       pars[0] = "mass";
       pars[1] = "mu1";
       // ----------------------------------------------//
  
       rooContainer->AddGenericPdf("exp",
         "exp((@0)*(@1))",pars);

       // make the vector of parameters as dependants of
       // the function to be fitted
       // ----------------------------------------------//
       std::vector<const char*> pars2(2,"t");
       pars2[0] = "mass";
       pars2[1] = "mu2";
       // ----------------------------------------------//
  
       rooContainer->AddGenericPdf("exp2",
         "exp((@0)*(@1))",pars2);
       // ----------------------------------------------//  
       std::vector<const char*> funs(2,"t");
       funs[0] = "exp";
       funs[1] = "exp2";
       // ----------------------------------------------//  
       rooContainer->ComposePdf("2xPdf","e1+e2",funs);

       // make a pdf for signal, (just an exp for now)
       std::vector<const char*> pars3(2,"t");
       pars3[0] = "mass_sig";
       pars3[1] = "mu3";

       rooContainer->AddGenericPdf("exp_sig",
	 "exp((@0)*(@1))",pars3);

       // ----------------------------------------------//

       rooContainer->CreateDataSet("mass");
       rooContainer->CreateDataSet("mass_sig");
    }
}

void LoopAll::TermRealPhotonAnalysis(int typerun) 
{
   if (typerun == 3){
 
     rooContainer->FitToData("2xPdf","mass",100,115,125,200);
     rooContainer->FitToData("exp_sig","mass_sig");
 
   }

}


void LoopAll::myFillHistPhotonAnalysis(Util *ut, int jentry) {

  if(PADEBUG) 
    cout << "myFillHist START"<<endl;

  counters[0]++;
  int histVal = ut->type2HistVal[ut->datatype[ut->current]];

  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    histoContainer[histVal]->Fill("pho_pt", p4->Pt());
  }
  
  Int_t in_endcap = 0;
  Float_t best_mass = 0;
  for (int i=0; i<pho_n-1; i++) {
    TLorentzVector *pg1= (TLorentzVector *) pho_p4->At(i);
    if (fabs(pg1->Eta()) > 1.479)
      in_endcap = 1;

    for (int j=i+1; j<pho_n; j++) {
      TLorentzVector *pg2= (TLorentzVector *) pho_p4->At(j);
      if (fabs(pg2->Eta()) > 1.479)
	in_endcap = 1;
      TLorentzVector higgs = (*pg1) + (*pg2);
      Float_t mass = higgs.M();
      if (mass > best_mass)
	best_mass = mass;
    }
  }     

  if (best_mass != 0) 
    histoContainer[histVal]->Fill("invmass", best_mass);
  
  if(PADEBUG) 
    cout<<"myFillHist END"<<endl;
}


void LoopAll::myFillHistPhotonAnalysisRed(Util * ut, int jentry) {


  if(PADEBUG) 
    cout << "myFillHistRed START"<<endl;

  int histVal = ut->type2HistVal[ut->datatype[ut->current]];
  counters[0]++;

  if (PADEBUG) 
    cout << "Getting Entries for Branches" << endl;
// From Here is the Standard Dec Review Selection/ gen Level studies
// and double sid-band distribution histograms
// Please do not edit!
  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_calopos->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_r9->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumptsolidconedr03->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_trksumpthollowconedr03->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  b_gen_n->GetEntry(jentry);
  b_gen_p4->GetEntry(jentry);
  b_gen_pdgid->GetEntry(jentry);
  b_gen_status->GetEntry(jentry);

  std::vector<Elec> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;
  for (int i=0; i<pho_n; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());
    //PreSelection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.5
       && pho_trksumpthollowconedr03[i] < 2*(3.5 + 0.001*pt)
       && pho_ecalsumetconedr03[i] < 2*(4.2 + 0.006*pt)
       && pho_hcalsumetconedr03[i] < 2*(2.2 + 0.0025*pt)
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
         candidate.r9 		= pho_r9[i];
         preselected_photons.push_back(candidate);
       }
  }

  if (PADEBUG) 
    cout << "Preselected Candidates" << endl;
  //Event Selection
  int n_preselected_pho = preselected_photons.size();
  //histoContainer[histVal]->Fill("h_n_sel",n_preselected_pho);

  // Sort Photons into Pt Order
  std::sort(preselected_photons.begin()
           ,preselected_photons.end()
           ,ElecP4greater); 
 
// Regular Event Selection begins here
  float best_mass = 0.;
  float best_pt = -1;
  int category=1;
  float min_r9  =1.;
  float max_eta =1.;

  if (PADEBUG) 
    cout << "Filling N-1 Plots" << endl;
  if (n_preselected_pho > 1 ){

  // Now The N-1 Plots --------------------------------------------------------------------
   std::vector<Elec>::iterator it_pho;

   for (it_pho  = preselected_photons.begin()
      ;it_pho != preselected_photons.end() 
      ;it_pho++ ) {
    
     bool pass_hoe    	= false;
     bool pass_sieie  	= false;
     bool pass_trkiso 	= false;
     bool pass_ecaliso  = false;
     bool pass_hcaliso  = false;
     bool is_eb 	= false;
     bool is_ee		= false;

     if (it_pho->hoe < 0.02) pass_hoe = true;
     if (((it_pho->sieie < 0.01)  && (fabs(it_pho->calopos->Eta()) < 1.4442)) 
         || (( it_pho->sieie < 0.028)
         && ((fabs(it_pho->calopos->Eta()) < 2.5) 
         && (fabs(it_pho->calopos->Eta()) > 1.566))) )      pass_sieie = true;
     if (it_pho->trkIso < (1.5 + 0.001*it_pho->p4->Pt()))   pass_trkiso = true;
     if (it_pho->ecalIso < (2.0 + 0.006*it_pho->p4->Pt()))  pass_ecaliso = true;
     if (it_pho->hcalIso < (2.0 + 0.0025*it_pho->p4->Pt())) pass_hcaliso = true;

     is_eb = (fabs(it_pho->calopos->Eta()) < 1.4442);
     is_ee = (fabs(it_pho->calopos->Eta()) < 2.5 && fabs(it_pho->calopos->Eta())>1.566);

  	 
     if (pass_sieie && pass_trkiso && pass_ecaliso && pass_hcaliso)
	histoContainer[histVal]->Fill("h_pho_hoe",it_pho->hoe);
     if (pass_hoe && pass_trkiso && pass_ecaliso && pass_hcaliso && is_eb)
	histoContainer[histVal]->Fill("h_pho_sieie_eb",it_pho->sieie);
     if (pass_hoe && pass_trkiso && pass_ecaliso && pass_hcaliso && is_ee)
	histoContainer[histVal]->Fill("h_pho_sieie_ee",it_pho->sieie);
     if (pass_sieie && pass_hoe && pass_ecaliso && pass_hcaliso);
	histoContainer[histVal]->Fill("h_pho_trksumpthollowconedr04"
				     ,it_pho->trkIso - 0.001*(it_pho->p4->Pt()));
     if (pass_sieie && pass_trkiso && pass_hoe && pass_hcaliso)
	histoContainer[histVal]->Fill("h_pho_ecalsumetconedr04"
				     ,it_pho->ecalIso - 0.006*(it_pho->p4->Pt()));
     if (pass_sieie && pass_trkiso && pass_hoe && pass_ecaliso)
	histoContainer[histVal]->Fill("h_pho_hcalsumetconedr04"
				     ,it_pho->hcalIso - 0.00025*(it_pho->p4->Pt()));

     if (pass_hoe && pass_sieie && pass_trkiso && pass_ecaliso && pass_hcaliso){
	histoContainer[histVal]->Fill("h_pho_pt",it_pho->p4->Pt());
	histoContainer[histVal]->Fill("h_pho_eta",fabs(it_pho->calopos->Eta()));	

     }
   
    }
 // -------------------------------------------------------------------------------------------
     bool dr_lead_match=false;
     bool dr_nlead_match=false;

     int closest_reco=-1;
     Elec leading, nleading;


     leading  = preselected_photons[0];
     nleading = preselected_photons[1];


  	 if (leading.p4->Pt() > 40.){
         
         // Dr gen-photon for leading and subleading candidate
         TLorentzVector *g_p4;
         for (int k=0;k<gen_n;k++){
           if (gen_pdgid[k] == 22 && gen_status[k] == 1 ){
	     g_p4 = (TLorentzVector*) gen_p4->At(k);
	     float dr1 =DeltaR(leading.p4->Phi(),g_p4->Phi()
	       	       	     ,leading.p4->Eta(),g_p4->Eta());
	     if (dr1 < 0.03) dr_lead_match = true;

	     float dr2 =DeltaR(nleading.p4->Phi(),g_p4->Phi()
		       	     ,nleading.p4->Eta(),g_p4->Eta());
	     if (dr2 < 0.03) dr_nlead_match = true;

	     if (dr1 < 0.05) histoContainer[histVal]->Fill("h_dr_lead",dr1);
	     if (dr2 < 0.05) histoContainer[histVal]->Fill("h_dr_nlead",dr2);
		
             if (dr1 < dr2 && dr_lead_match) closest_reco = 1;
	     else if (dr2 < dr1 && dr_nlead_match) closest_reco = 2;
	     
 	     if (dr_lead_match || dr_nlead_match) break;
	    }
         }
	 

         if (PADEBUG) 
    	   cout << "Looking for a Higgs" << endl;
         // -------------------------------------------------------
         TLorentzVector Higgs = (*(preselected_photons[0].p4))
                                 +(*(preselected_photons[1].p4));
           float mass = Higgs.M();
           float h_pt = Higgs.Pt();
           // Determine the Category of the event
           // -> Which histogram is filled
           min_r9  = min(leading.r9,nleading.r9);
	   max_eta = max(fabs(leading.calopos->Eta())
		      ,fabs(nleading.calopos->Eta()));
	   if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
	   if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
	   if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
	   if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;

           if (mass > 100. && mass < 150.){
             //Good event, passes preselection and acceptance cuts

             int pass_selection[2];
             int pass_isolation[2];
             int in_iso_gap[2];
    
            //Now do selection on leading photon
             pass_selection[0] = leading.hoe < 0.02
                	       && (((leading.sieie < 0.01)  && (fabs(leading.calopos->Eta()) < 1.4442)) 
                		  || (( leading.sieie < 0.028)
                   	       && ((fabs(leading.calopos->Eta()) < 2.5) && (fabs(leading.calopos->Eta()) > 1.566))) );
             pass_isolation[0] =  	
                                 leading.trkIso < (1.5 + 0.001*leading.p4->Pt())	
                               && leading.hcalIso < (2.0 + 0.0025*leading.p4->Pt())
                	       && leading.ecalIso < (2.0 + 0.006*leading.p4->Pt());
	     in_iso_gap[0]     = //false;	
				 !pass_isolation[0]
			        && leading.trkIso < (3 + 0.001*leading.p4->Pt());
	 
             //Selection on next to leading photon
             pass_selection[1] = nleading.hoe < 0.02
                	       && (((nleading.sieie < 0.01)  && (fabs(nleading.calopos->Eta()) < 1.4442)) 
                		  || (( nleading.sieie < 0.028)
                  	       && ((fabs(nleading.calopos->Eta()) < 2.5) && (fabs(nleading.calopos->Eta()) > 1.566))) );
             pass_isolation[1] = 
                	          nleading.trkIso < (1.5 + 0.001*nleading.p4->Pt())
                               && (nleading.hcalIso < (2.0 + 0.0025*nleading.p4->Pt()))
                	       && (nleading.ecalIso < (2.0 + 0.006*nleading.p4->Pt()));
	     in_iso_gap[1]     = //false; 
			       !pass_isolation[1]
			       && nleading.trkIso < (3 + 0.001*nleading.p4->Pt());

             //Correlation Plots
	     /*
             if ( closest_reco == 2 ){
               histoContainer[histVal]->Fill("h_corr_trckiso_vs_hoe",leading.trkIso - 0.001*leading.p4->Pt()
								  ,leading.hoe);
               histoContainer[histVal]->Fill("h_corr_trckiso_vs_sieie",leading.trkIso - 0.001*leading.p4->Pt()
								  ,leading.sieie);

               histoContainer[histVal]->Fill("h_corr_ecaliso_vs_hoe",leading.ecalIso - 0.006*leading.p4->Pt()
								  ,leading.hoe);
               histoContainer[histVal]->Fill("h_corr_ecaliso_vs_sieie",leading.ecalIso - 0.006*leading.p4->Pt()
								  ,leading.sieie);

               histoContainer[histVal]->Fill("h_corr_hcaliso_vs_hoe",leading.hcalIso - 0.0025*leading.p4->Pt()
								  ,leading.hoe);
               histoContainer[histVal]->Fill("h_corr_hcaliso_vs_sieie",leading.hcalIso - 0.0025*leading.p4->Pt()
	      						  	  ,leading.sieie);
	     }  
             if ( closest_reco ==1){
               histoContainer[histVal]->Fill("h_corr_trckiso_vs_hoe",nleading.trkIso - 0.001*nleading.p4->Pt()
								  ,nleading.hoe);
               histoContainer[histVal]->Fill("h_corr_trckiso_vs_sieie",nleading.trkIso - 0.001*nleading.p4->Pt()
								  ,nleading.sieie);

               histoContainer[histVal]->Fill("h_corr_ecaliso_vs_hoe",nleading.ecalIso - 0.006*nleading.p4->Pt()
								  ,nleading.hoe);
               histoContainer[histVal]->Fill("h_corr_ecaliso_vs_sieie",nleading.ecalIso - 0.006*nleading.p4->Pt()
								  ,nleading.sieie);

               histoContainer[histVal]->Fill("h_corr_hcaliso_vs_hoe",nleading.hcalIso - 0.0025*nleading.p4->Pt()
								  ,nleading.hoe);
               histoContainer[histVal]->Fill("h_corr_hcaliso_vs_sieie",nleading.hcalIso - 0.0025*nleading.p4->Pt()
	      						  	  ,nleading.sieie);
	     }  
	     //Double Sideband Method

	*/
	     if (pass_selection[0])  histoContainer[histVal]->Fill(
					"h_pho_trksumpthollowconedr04_leading_pass"
					,leading.trkIso - 0.001*leading.p4->Pt());

	     if (!pass_selection[0]) histoContainer[histVal]->Fill(
					"h_pho_trksumpthollowconedr04_leading_fail"
					,leading.trkIso - 0.001*leading.p4->Pt());
	     if (pass_selection[0] && pass_isolation[0]){

	       if (pass_selection[1])  histoContainer[histVal]->Fill(
					"h_pho_trksumpthollowconedr04_nleading_pass"
					,nleading.trkIso - 0.001*nleading.p4->Pt());

	       if (!pass_selection[1]) histoContainer[histVal]->Fill(
					"h_pho_trksumpthollowconedr04_nleading_fail"
					,nleading.trkIso - 0.001*nleading.p4->Pt());
	     }

	     // Fill Side-bands and Signal Regions
	     if (!in_iso_gap[0]){
              histoContainer[histVal]->Fill("h_sideband_leading",
                                         pass_isolation[0],pass_selection[0]);
              if (category == 1)      histoContainer[histVal]->Fill("h_sideband_leading_cat1",
                                         pass_isolation[0],pass_selection[0]);
              else if (category == 2) histoContainer[histVal]->Fill("h_sideband_leading_cat2",
                                         pass_isolation[0],pass_selection[0]);
              else if (category == 3) histoContainer[histVal]->Fill("h_sideband_leading_cat3",
                                         pass_isolation[0],pass_selection[0]);
              else if (category == 4) histoContainer[histVal]->Fill("h_sideband_leading_cat4",
                                         pass_isolation[0],pass_selection[0]);
              
	      if (pass_selection[0] && pass_isolation[0] && !in_iso_gap[1]){
               histoContainer[histVal]->Fill("h_sideband_nleading",
                                        pass_isolation[1],pass_selection[1]);
               if (category == 1)     histoContainer[histVal]->Fill("h_sideband_nleading_cat1",
                                        pass_isolation[1],pass_selection[1]);
               else if (category == 2)histoContainer[histVal]->Fill("h_sideband_nleading_cat2",
                                        pass_isolation[1],pass_selection[1]);
               else if (category == 3)histoContainer[histVal]->Fill("h_sideband_nleading_cat3",
                                        pass_isolation[1],pass_selection[1]);
               else if (category == 4)histoContainer[histVal]->Fill("h_sideband_nleading_cat4",
                                        pass_isolation[1],pass_selection[1]);

               if (pass_selection[1] && pass_isolation[1]){
		 histoContainer[histVal]->Fill("h_pho_pt",leading.p4->Pt());
		 histoContainer[histVal]->Fill("h_pho_pt",nleading.p4->Pt());
                 best_mass = mass;
 		 best_pt   = h_pt;

                 if (PADEBUG) 
    	          cout << "Higgs Candidate Found in Event" << endl;
               }
              }
	    }

	     // Histogram used to calculate A = Njy/Nyj
	     if(   pass_selection[0]
		&& pass_selection[1]
		&& pass_isolation[0]
		&& pass_isolation[1]
		&& closest_reco > 0){
		 if (closest_reco == 1){
		     histoContainer[histVal]->Fill("h_selected_asymm",0.);
		  if(category ==1) histoContainer[histVal]->Fill("h_selected_asymm_cat1",0.);
		  if(category ==2) histoContainer[histVal]->Fill("h_selected_asymm_cat2",0.);
		  if(category ==3) histoContainer[histVal]->Fill("h_selected_asymm_cat3",0.);
		  if(category ==4) histoContainer[histVal]->Fill("h_selected_asymm_cat4",0.);
		 }
		 if (closest_reco ==2){
		     histoContainer[histVal]->Fill("h_selected_asymm",1.);
		  if(category ==1) histoContainer[histVal]->Fill("h_selected_asymm_cat1",1.);
		  if(category ==2) histoContainer[histVal]->Fill("h_selected_asymm_cat2",1.);
		  if(category ==3) histoContainer[histVal]->Fill("h_selected_asymm_cat3",1.);
		  if(category ==4) histoContainer[histVal]->Fill("h_selected_asymm_cat4",1.);
		 }
		}
           }
     }
   }

  if( best_mass > 90 && best_mass < 210){ 
    histoContainer[histVal]->Fill("h_mass", best_mass);
    if(category ==1) histoContainer[histVal]->Fill("h_mass_cat1", best_mass);
    if(category ==2) histoContainer[histVal]->Fill("h_mass_cat2", best_mass);
    if(category ==3) histoContainer[histVal]->Fill("h_mass_cat3", best_mass);
    if(category ==4) histoContainer[histVal]->Fill("h_mass_cat4", best_mass);
  }

  histoContainer[histVal]->Fill("h_pt", best_pt);
  if(category ==1) histoContainer[histVal]->Fill("h_pt_cat1", best_pt);
  if(category ==2) histoContainer[histVal]->Fill("h_pt_cat2", best_pt);
  if(category ==3) histoContainer[histVal]->Fill("h_pt_cat3", best_pt);
  if(category ==4) histoContainer[histVal]->Fill("h_pt_cat4", best_pt);

  
  if(PADEBUG) 
    cout<<"myFillHistRed END"<<endl;
}

void LoopAll::myReducePhotonAnalysis(Util * ut, int jentry) {

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis START"<<endl;

  //count all events
  countersred[0]++;

  if(outputFile) {
    if(makeOutputTree) {
      
      //first selection and fill output tree
      if(!myFillReducedVarPhotonAnalysis(ut, jentry)) 
	return;
      
      //additional selection
      if(!mySelectEventRedPhotonAnalysis(ut, jentry)) 
	return;

      countersred[1]++;

      outputEvents++;
      if(PADEBUG) 
	cout<<"before fill"<<endl;

      outputTree->Fill();
      if(PADEBUG) 
	cout<<"after fill"<<endl;

      if(outputEvents==100) {
	outputEvents=0;
	outputTree->Write(0,TObject::kWriteDelete);
      }
    }
  }

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis END"<<endl;
}

void LoopAll::myGetBranchPhotonAnalysis() {
  b_pho_n = fChain->GetBranch("pho_n");
  b_pho_p4 = fChain->GetBranch("pho_p4");
  b_pho_calopos = fChain->GetBranch("pho_calopos");
  b_pho_hoe = fChain->GetBranch("pho_hoe");
  b_pho_r9 = fChain->GetBranch("pho_r9");
  b_pho_sieie = fChain->GetBranch("pho_sieie");
  b_pho_ecalsumetconedr03 = fChain->GetBranch("pho_ecalsumetconedr03");
  b_pho_ecalsumetconedr04 = fChain->GetBranch("pho_ecalsumetconedr04");
  b_pho_hcalsumetconedr03 = fChain->GetBranch("pho_hcalsumetconedr03");
  b_pho_hcalsumetconedr04 = fChain->GetBranch("pho_hcalsumetconedr04");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumpthollowconedr04 = fChain->GetBranch("pho_trksumpthollowconedr04");
  b_pho_trksumpthollowconedr03 = fChain->GetBranch("pho_trksumpthollowconedr03");
  b_pho_isEB = fChain->GetBranch("pho_isEB");
  b_pho_isEE = fChain->GetBranch("pho_isEE");
  b_pho_haspixseed = fChain->GetBranch("pho_haspixseed");

  b_gen_n = fChain->GetBranch("gp_n");
  b_gen_p4 = fChain->GetBranch("gp_p4");
  b_gen_status = fChain->GetBranch("gp_status");
  b_gen_pdgid = fChain->GetBranch("gp_pdgid");

}

int LoopAll::myFillReducedVarPhotonAnalysis(Util * ut, int jentry) {
  if(PADEBUG) 
    cout<<"myFillReduceVar START"<<endl;
  
  b_pho_p4->GetEntry(jentry);
  for (int i=0; i<pho_n-1; i++) {
    TLorentzVector *pho_temp = (TLorentzVector *) pho_p4->At(i);
    pho_Et[i]=pho_temp->Et();
  }

  if(PADEBUG) 
    cout<<"myFillReduceVar END"<<endl;

}
//This relates to branchs being read in
void LoopAll::mySetBranchAddressRedPhotonAnalysis() {
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
  fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
  fChain->SetBranchAddress("pho_hoe", &pho_hoe, &b_pho_hoe);
  fChain->SetBranchAddress("pho_r9", &pho_r9, &b_pho_r9);
  fChain->SetBranchAddress("pho_sieie", &pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
  fChain->SetBranchAddress("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
  fChain->SetBranchAddress("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
  fChain->SetBranchAddress("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
  fChain->SetBranchAddress("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
  fChain->SetBranchAddress("pho_trksumptsolidconedr04", &pho_trksumptsolidconedr04, &b_pho_trksumptsolidconedr04);
  fChain->SetBranchAddress("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  fChain->SetBranchAddress("pho_trksumpthollowconedr03", &pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
  fChain->SetBranchAddress("pho_isEB", &pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", &pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_haspixseed", &pho_haspixseed, &b_pho_haspixseed);

  fChain->SetBranchAddress("gp_n", &gen_n, &b_gen_n);
  fChain->SetBranchAddress("gp_p4", &gen_p4, &b_gen_p4);
  fChain->SetBranchAddress("gp_status", &gen_status, &b_gen_status);
  fChain->SetBranchAddress("gp_pdgid", &gen_pdgid, &b_gen_pdgid);
}


int LoopAll::mySelectEventRedPhotonAnalysis(Util * ut, int jentry) {
  
  // preselection at the end
  int selectevent=0;

  b_pho_n->GetEntry(jentry);


  if (pho_n > 1)
    selectevent = 1;
  else
    selectevent = 0;

  return selectevent;
}
