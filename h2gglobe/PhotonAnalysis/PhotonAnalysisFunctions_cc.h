#define PADEBUG 0


void LoopAll::InitRealPhotonAnalysis(Util *ut, int typerun) {

  // Book histos only if not reduce step
  if (typerun != 1) {    
   for (int ind=0;ind<ut->ntypes;ind++){
    histoContainer[ind]->Add("pho_pt", 100, 0, 100);
    histoContainer[ind]->Add("invmass_barrel", 200, 0, 200);
    histoContainer[ind]->Add("invmass_endcap", 200, 0, 200);
    histoContainer[ind]->Add("h_mass","M #gamma #gamma - GeV/c^{2}", 120, 80, 200);
    histoContainer[ind]->Add("h_pho_hoe", 100, 0, 0.1);
    histoContainer[ind]->Add("h_pho_sieie", 100, 0, 0.07);
    histoContainer[ind]->Add("h_pho_ecalsumetconedr03", 100, 0, 20);
    histoContainer[ind]->Add("h_pho_ecalsumetconedr04", 100, 0, 20);
    histoContainer[ind]->Add("h_pho_hcalsumetconedr03", 100, 0, 20);
    histoContainer[ind]->Add("h_pho_hcalsumetconedr04", 100, 0, 20);
    histoContainer[ind]->Add("h_pho_trksumptsolidconedr03", 100, 0, 20);
    histoContainer[ind]->Add("h_pho_trksumptsolidconedr04", 100, 0, 20);
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

  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumptsolidconedr03->GetEntry(jentry); 
  b_pho_trksumptsolidconedr04->GetEntry(jentry);  
  b_pho_trksumpthollowconedr03->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry); 

  std::vector<TLorentzVector *> selected_photons;  

  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    histoContainer[histVal]->Fill("pho_pt", p4->Pt());
    histoContainer[histVal]->Fill("h_pho_hoe", pho_hoe[i]);
    histoContainer[histVal]->Fill("h_pho_sieie", pho_sieie[i]);
    histoContainer[histVal]->Fill("h_pho_ecalsumetconedr03", pho_ecalsumetconedr03[i]);
    histoContainer[histVal]->Fill("h_pho_ecalsumetconedr04", pho_ecalsumetconedr04[i]);
    histoContainer[histVal]->Fill("h_pho_hcalsumetconedr03", pho_hcalsumetconedr03[i]);
    histoContainer[histVal]->Fill("h_pho_hcalsumetconedr04", pho_hcalsumetconedr04[i]);
    histoContainer[histVal]->Fill("h_pho_trksumptsolidconedr03", pho_trksumptsolidconedr03[i]);
    histoContainer[histVal]->Fill("h_pho_trksumptsolidconedr04", pho_trksumptsolidconedr04[i]);
  
    float pt  = p4->Pt();
    float eta = fabs(p4->Eta());

    if(
	(! pho_haspixseed[i])
     && pt > 30. 
     && pho_hoe[i] < 0.02
     && pho_trksumpthollowconedr04[i] < (1.5 + 0.001*pt)
     && pho_ecalsumetconedr04[i] < (2.0 + 0.006*pt)
     && pho_hcalsumetconedr04[i] < (2.0 + 0.0025*pt)
     && (   ((pho_sieie[i] < 0.01) && (eta < 1.44)) 
	 || ((pho_sieie[i] <0.028) && ((eta < 2.5) && (eta > 1.57))) 
	)
     && (eta < 1.44) || ((eta > 1.57) && (eta < 2.5)) 
 
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
  b_pho_hoe = fChain->GetBranch("pho_hoe");
  b_pho_sieie = fChain->GetBranch("pho_sieie");
  b_pho_ecalsumetconedr03 = fChain->GetBranch("pho_ecalsumetconedr03");
  b_pho_ecalsumetconedr04 = fChain->GetBranch("pho_ecalsumetconedr04");
  b_pho_hcalsumetconedr03 = fChain->GetBranch("pho_hcalsumetconedr03");
  b_pho_hcalsumetconedr04 = fChain->GetBranch("pho_hcalsumetconedr04");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  //b_pho_trksumpthollowconedr04 = fChain->GetBranch("pho_trksumpthollowconedr04");
  //b_pho_trksumpthollowconedr04 = fChain->GetBranch("pho_trksumpthollowconedr04");
  b_pho_isEB = fChain->GetBranch("pho_isEB");
  b_pho_isEE = fChain->GetBranch("pho_isEE");
  b_pho_haspixseed = fChain->GetBranch("pho_haspixseed");
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
  fChain->SetBranchAddress("pho_hoe", &pho_hoe, &b_pho_hoe);
  fChain->SetBranchAddress("pho_sieie", &pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
  fChain->SetBranchAddress("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
  fChain->SetBranchAddress("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
  fChain->SetBranchAddress("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
  fChain->SetBranchAddress("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
  fChain->SetBranchAddress("pho_trksumptsolidconedr04", &pho_trksumptsolidconedr04, &b_pho_trksumptsolidconedr04);
  //fChain->SetBranchAddress("pho_trksumpthollowconedr03", &pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
  //fChain->SetBranchAddress("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  fChain->SetBranchAddress("pho_isEB", &pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", &pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_haspixseed", &pho_haspixseed, &b_pho_haspixseed);
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

