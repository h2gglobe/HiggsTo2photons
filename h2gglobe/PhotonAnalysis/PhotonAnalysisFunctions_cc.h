#define PADEBUG 0


void LoopAll::InitRealPhotonAnalysis(Util *ut, int typerun) {

  // Book histos only if not reduce step
  if (typerun != 1) {    
   for (int ind=0;ind<ut->ntypes;ind++){
//    histoContainer[ind]->Add("pho_pt", 100, 0, 100);
//    histoContainer[ind]->Add("invmass_barrel", 200, 0, 200);
//    histoContainer[ind]->Add("invmass_endcap", 200, 0, 200);
    histoContainer[ind]->Add("h_mass","M #gamma #gamma - GeV/c^{2}", 10, 100, 150);
//    histoContainer[ind]->Add("h_pho_hoe", 100, 0, 0.1);
//    histoContainer[ind]->Add("h_pho_sieie", 100, 0, 0.07);
//    histoContainer[ind]->Add("h_pho_ecalsumetconedr03", 100, 0, 20);
//    histoContainer[ind]->Add("h_pho_ecalsumetconedr04", 100, 0, 20);
//    histoContainer[ind]->Add("h_pho_hcalsumetconedr03", 100, 0, 20);
//    histoContainer[ind]->Add("h_pho_hcalsumetconedr04", 100, 0, 20);
//    histoContainer[ind]->Add("h_pho_trksumptsolidconedr03", 100, 0, 20);
//    histoContainer[ind]->Add("h_pho_trksumpthollowconedr04", 100, 0, 20);
    histoContainer[ind]->Add("h_sideband_leading",  2, -0.5, 1.5, 2, -0.5, 1.5);
    histoContainer[ind]->Add("h_sideband_nleading", 2, -0.5, 1.5, 2, -0.5, 1.5);
   }
  }

    if (typerun == 3){
 
       rooContainer->AddRealVar("mass",100.,200.);
       rooContainer->AddRealVar("mu1",-0.04,-1.0,-0.001);
       rooContainer->AddRealVar("mu2",-0.04,-1.0,-0.001);
 
       // make the vector of parameters as dependants of
       // the function to be fitted
       // ----------------------------------------------//
       std::vector<const char*> pars(3,"t");
       pars[0] = "mass";
       pars[1] = "mu1";
       pars[2] = "mu2";
       // ----------------------------------------------//
  
       rooContainer->AddGenericPdf("exp",
         "exp((@0)*(@1))+exp((@0)*(@2))",pars);
  
       rooContainer->CreateDataSet("mass");
    }
}

void LoopAll::TermRealPhotonAnalysis(int typerun) 
{
   if (typerun == 3){
 
     rooContainer->FitToData("exp","mass");
 
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
 //   histoContainer[histVal]->Fill("pho_pt", p4->Pt());
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
//  b_pho_trksumptsolidconedr04->GetEntry(jentry);  
 // b_pho_trksumpthollowconedr03->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  b_gp_n->GetEntry();
  b_gp_p4->GetEntry();
  b_gp_status->GetEntry();
  b_gp_pdgid->GetEntry();

//  struct Elec{
//    TLorentzVector *p4;
//    bool pixSeed;
//    double trkIso;
//    double ecalIso;
//    double hcalIso;
//    double sieie;
//    double hoe;
//  };

  std::vector<Elec> preselected_photons;  
  std::vector<Elec> selected_photons; 
 
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    TLorentzVector *calopos = (TLorentzVector *) pho_calopos->At(i);
//    histoContainer[histVal]->Fill("pho_pt", p4->Pt());
//    histoContainer[histVal]->Fill("h_pho_hoe", pho_hoe[i]);
//    histoContainer[histVal]->Fill("h_pho_sieie", pho_sieie[i]);
//    histoContainer[histVal]->Fill("h_pho_ecalsumetconedr03", pho_ecalsumetconedr03[i]);
//    histoContainer[histVal]->Fill("h_pho_ecalsumetconedr04", pho_ecalsumetconedr04[i]);
//    histoContainer[histVal]->Fill("h_pho_hcalsumetconedr03", pho_hcalsumetconedr03[i]);
//    histoContainer[histVal]->Fill("h_pho_hcalsumetconedr04", pho_hcalsumetconedr04[i]);
//    histoContainer[histVal]->Fill("h_pho_trksumptsolidconedr03", pho_trksumptsolidconedr03[i]);
//    histoContainer[histVal]->Fill("h_pho_trksumpthollowconedr04", pho_trksumpthollowconedr04[i]);
    float pt  = p4->Pt();
    float eta = fabs(p4->Eta());
    //PreSelection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.10
       //&& pho_trksumpthollowconedr04[i] < (7.0 + 0.002*pt)
       && pho_ecalsumetconedr04[i] < (8.4 + 0.012*pt)
       && pho_hcalsumetconedr04[i] < (4.4 + 0.005*pt)
       &&((eta < 1.4442) || ((eta > 1.566) && (eta < 2.5))) 
       ) {
         Elec candidate;
         candidate.p4 = p4;
	 candidate.calopos = calopos;
         candidate.pixSeed = pho_haspixseed[i];
         candidate.trkIso = pho_trksumpthollowconedr04[i];
         candidate.ecalIso = pho_ecalsumetconedr04[i];
         candidate.hcalIso = pho_hcalsumetconedr04[i];
         candidate.sieie = pho_sieie[i];
         candidate.hoe = pho_hoe[i];
         preselected_photons.push_back(candidate);
       }
  }

  float best_mass = 0.;
  //Event Selection
  int n_preselected_pho = preselected_photons.size();
  if (n_preselected_pho > 1 ){
    for (int i=0; i< n_preselected_pho-1; i++){
      for (int j=i+1;j<n_preselected_pho;j++){
        if (preselected_photons[i].p4->Pt() > 40. 
	   || preselected_photons[j].p4->Pt() > 40.){
          
           TLorentzVector Higgs = (*(preselected_photons[i].p4))
                                 +(*(preselected_photons[j].p4));
           float mass = Higgs.M();
             if (mass > 100. && mass < 150.){
             //Good event, passes preselection and acceptance cuts
             Elec leading  = preselected_photons[i];
             Elec nleading = preselected_photons[j];

             int pass_selection[2];
             int pass_isolation[2];
    
            //Now do selection on leading photon
             pass_selection[0] = (!leading.pixSeed)
                && leading.hoe < 0.02
                && leading.ecalIso < (2.0 + 0.006*leading.p4->Pt())
                && leading.hcalIso < (2.0 + 0.0025*leading.p4->Pt())
                && (((leading.sieie < 0.01)  && (leading.p4->Eta() < 1.4442)) 
                || (( leading.sieie < 0.028)
                   && ((leading.calopos->Eta() < 2.5) && (leading.calopos->Eta() > 1.566))) );
             pass_isolation[0] =  leading.trkIso < (1.5 + 0.001*leading.p4->Pt());
             //Selection on next to leading photon
             pass_selection[1] = (!nleading.pixSeed)
                && nleading.hoe < 0.02
                && nleading.ecalIso < (2.0 + 0.006*nleading.p4->Pt())
                && nleading.hcalIso < (2.0 + 0.0025*nleading.p4->Pt())
                && (((nleading.sieie < 0.01)  && (nleading.p4->Eta() < 1.4442)) 
                || (( nleading.sieie < 0.028)
                  && ((nleading.calopos->Eta() < 2.5) && (nleading.calopos->Eta() > 1.566))) );
             pass_isolation[1] =  nleading.trkIso < (1.5 + 0.001*nleading.p4->Pt());
             //Double Sideband Method
             histoContainer[histVal]->Fill("h_sideband_leading",
                                         pass_isolation[0],pass_selection[0]);
             if (pass_selection[0] && pass_isolation[0]){
               histoContainer[histVal]->Fill("h_sideband_nleading",
                                         pass_isolation[1],pass_selection[1]);
               if (pass_selection[1] && pass_isolation[1]){
                 best_mass = mass;
               }
             }
           }
         } 
       }
     }
   }

     
  histoContainer[histVal]->Fill("h_mass", best_mass);
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
  b_pho_calopos = fChain->GetBranch("pho_p4");
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
