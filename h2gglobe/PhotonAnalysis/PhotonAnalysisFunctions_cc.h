#define PADEBUG 0


void LoopAll::InitRealPhotonAnalysis(Util *ut, int typerun) {

  // Book histos only if not reduce step
  if (typerun != 1) {    
   for (int ind=0;ind<ut->ntypes;ind++){
    histoContainer[ind]->Add("h_mass","M #gamma #gamma - GeV/c^{2}", 10, 100, 150);
    histoContainer[ind]->Add("h_pt","p_{T} #gamma #gamma - GeV/c", 50, 0, 150);
    histoContainer[ind]->Add("h_n_sel","pre-selected photons", 4, 0, 4);
    histoContainer[ind]->Add("h_pho_hoe","H/E", 100, 0, 0.05);
    histoContainer[ind]->Add("h_pho_pt","p_{T} ^{#gamma}", 70, 30, 100);

    // Isolation Histograms for Pass/Fail distributions ---------------------------------------------------
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_leading_pass","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_leading_fail","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_nleading_pass","TrkIso hollow DR04", 20, 0, 3);
    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04_nleading_fail","TrkIso hollow DR04", 20, 0, 3);
    // ----------------------------------------------------------------------------------------------------

    histoContainer[ind]->Add("h_dr_lead",150,0,0.15);
    histoContainer[ind]->Add("h_dr_nlead",150,0,0.15);

    histoContainer[ind]->Add("h_sideband_leading",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading", 2, -0.5, 1.5, 2, -0.5, 1.5);

    histoContainer[ind]->Add("h_selected_asymm", 2,0,2);
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

// From Here is the Standard Dec Review Selection/ gen Level studies
// and double sid-band distribution histograms
// Please do not edit!
  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_calopos->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumptsolidconedr03->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  b_gp_n->GetEntry(jentry);
  b_gp_p4->GetEntry(jentry);
  b_gp_pdgid->GetEntry(jentry);
  b_gp_status->GetEntry(jentry);

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
  histoContainer[histVal]->Fill("h_n_sel",n_preselected_pho);

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


         // Dr gen-photon for leading and subleading candidate
         TLorentzVector *gen_p4;
         for (int k=0;k<gp_n;k++){
           if (gp_pdgid[k] == 22 && gp_status[k] == 1 ){
	     gen_p4 = (TLorentzVector*) gp_p4->At(k);
	     float dr1 =DeltaR(leading.p4->Phi(),gen_p4->Phi()
	       	       	     ,leading.p4->Eta(),gen_p4->Eta());
	     if (dr1 < 0.01) dr_lead_match = true;

	     float dr2 =DeltaR(nleading.p4->Phi(),gen_p4->Phi()
		       	     ,nleading.p4->Eta(),gen_p4->Eta());
	     if (dr2 < 0.01) dr_nlead_match = true;

	     if (dr1 < 0.15) histoContainer[histVal]->Fill("h_dr_lead",dr1);
	     if (dr2 < 0.15) histoContainer[histVal]->Fill("h_dr_nlead",dr2);
	    }
         }
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
             //Double Sideband Method

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
             histoContainer[histVal]->Fill("h_sideband_leading",
                                         pass_isolation[0],pass_selection[0]);
             if (pass_selection[0] && pass_isolation[0]){
               histoContainer[histVal]->Fill("h_sideband_nleading",
                                        pass_isolation[1],pass_selection[1]);
               if (pass_selection[1] && pass_isolation[1]){
		 histoContainer[histVal]->Fill("h_pho_pt",leading.p4->Pt());
		 histoContainer[histVal]->Fill("h_pho_pt",nleading.p4->Pt());
                 best_mass = mass;
 		 best_pt   = h_pt;
               }
             }

	     // Histogram used to calculate A = Njy/Nyj
	     if(   pass_selection[0]
		&& pass_selection[1]
		&& pass_isolation[0]
		&& pass_isolation[1]){
		 if (dr_lead_match && !dr_nlead_match)
		     histoContainer[histVal]->Fill("h_selected_asymm",0.);
		 if (dr_nlead_match && !dr_lead_match)
		     histoContainer[histVal]->Fill("h_selected_asymm",1.);
		}
           }
     }
   }

     
  histoContainer[histVal]->Fill("h_mass", best_mass);
  histoContainer[histVal]->Fill("h_pt", best_pt);


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
  b_pho_sieie = fChain->GetBranch("pho_sieie");
  b_pho_ecalsumetconedr03 = fChain->GetBranch("pho_ecalsumetconedr03");
  b_pho_ecalsumetconedr04 = fChain->GetBranch("pho_ecalsumetconedr04");
  b_pho_hcalsumetconedr03 = fChain->GetBranch("pho_hcalsumetconedr03");
  b_pho_hcalsumetconedr04 = fChain->GetBranch("pho_hcalsumetconedr04");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumpthollowconedr04 = fChain->GetBranch("pho_trksumpthollowconedr04");
  b_pho_isEB = fChain->GetBranch("pho_isEB");
  b_pho_isEE = fChain->GetBranch("pho_isEE");
  b_pho_haspixseed = fChain->GetBranch("pho_haspixseed");

  b_gp_n = fChain->GetBranch("gp_n");
  b_gp_p4 = fChain->GetBranch("gp_p4");
  b_gp_status = fChain->GetBranch("gp_status");
  b_gp_pdgid = fChain->GetBranch("gp_pdgid");

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
  fChain->SetBranchAddress("pho_sieie", &pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
  fChain->SetBranchAddress("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
  fChain->SetBranchAddress("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
  fChain->SetBranchAddress("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
  fChain->SetBranchAddress("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
  fChain->SetBranchAddress("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  //fChain->SetBranchAddress("pho_trksumpthollowconedr03", &pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
  //fChain->SetBranchAddress("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  fChain->SetBranchAddress("pho_isEB", &pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", &pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_haspixseed", &pho_haspixseed, &b_pho_haspixseed);

  fChain->SetBranchAddress("gp_n", &gp_n, &b_gp_n);
  fChain->SetBranchAddress("gp_p4", &gp_p4, &b_gp_p4);
  fChain->SetBranchAddress("gp_status", &gp_status, &b_gp_status);
  fChain->SetBranchAddress("gp_pdgid", &gp_pdgid, &b_gp_pdgid);
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
