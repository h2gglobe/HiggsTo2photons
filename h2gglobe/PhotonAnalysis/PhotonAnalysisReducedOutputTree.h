void LoopAll::PhotonAnalysisReducedOutputTree() {
  UtilInstance->outputTree->Branch("pho_n", &pho_n, "pho_n/I");
  UtilInstance->outputTree->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);
  UtilInstance->outputTree->Branch("pho_calopos", "TClonesArray", &pho_calopos, 32000, 0);
  UtilInstance->outputTree->Branch("pho_hoe",  &pho_hoe, "pho_hoe[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_sieie", &pho_sieie,"pho_sieie[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, "pho_ecalsumetconedr03[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, "pho_ecalsumetconedr04[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, "pho_hcalsumetconedr03[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, "pho_hcalsumetconedr04[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, "pho_trksumptsolidconedr03[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, "pho_trksumptsolidconedr03[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, "pho_trksumpthollowconedr04[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, "pho_trksumpthollowconedr04[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_isEB", &pho_isEB, "pho_isEB[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_isEE", &pho_isEE, "pho_isEE[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_haspixseed", &pho_haspixseed, "pho_haspixseed[pho_n]/F");
  UtilInstance->outputTree->Branch("pho_Et", &pho_Et, "pho_Et[pho_n]/F");
  
  //Generator Information
  
  UtilInstance->outputTree->Branch("gp_n",  &gp_n, "gp_n/I");
  UtilInstance->outputTree->Branch("gp_p4", "TClonesArray", &gp_p4,32000,0);
  UtilInstance->outputTree->Branch("gp_vtx", "TClonesArray", &gp_vtx,32000,0);
  UtilInstance->outputTree->Branch("gp_status",&gp_status, "gp_status/I");
  UtilInstance->outputTree->Branch("gp_pdgid",&gp_pdgid, "gp_pdgid/I");
  UtilInstance->outputTree->Branch("gp_mother",&gp_mother, "gp_mother/I");

}

