#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

GlobeLeptons::GlobeLeptons() {
  leptons_cut_eta_rec=15.;
  leptons_cut_et_rec=0.;
}

void GlobeLeptons::defineBranch(GlobeAnalyzer* ana) {
  lpt_p4 = new TClonesArray("TLorentzVector", MAX_LEPTONS);

  ana->Branch("lpt_n", &lpt_n, "lpt_n/I");
  ana->Branch("lpt_emu_n", &lpt_emu_n, "lpt_emu_n/I");
  ana->Branch("lpt_mu_n", &lpt_mu_n, "lpt_mu_n/I");
  ana->Branch("lpt_el_n", &lpt_el_n, "lpt_el_n/I");
  ana->Branch("lpt_pho_n", &lpt_pho_n, "lpt_pho_n/I");
  ana->Branch("lpt_pdgid", &lpt_pdgid, "lpt_pdgid[lpt_n]/I");
  ana->Branch("lpt_ind", &lpt_ind, "lpt_ind[lpt_n]/I");
  ana->Branch("lpt_duplicate", &lpt_duplicate, "lpt_duplicate[lpt_n]/I");
  ana->Branch("lpt_p4", "TClonesArray", &lpt_p4, 32000, 0);
  ana->Branch("lpt_indgen", &lpt_indgen, "lpt_indgen[lpt_n]/I");// filled in GlobeReducedGen
  ana->Branch("lpt_drmatch", &lpt_drmatch, "lpt_drmatch[lpt_n]/F");// filled in GlobeReducedGen
}

void GlobeLeptons::Zero() {
  lpt_n=0;
  lpt_mu_n=0;
  lpt_el_n=0;
  lpt_emu_n=0;
  lpt_pho_n=0;
  lpt_p4->Clear();  
}

void GlobeLeptons::addMuons(GlobeMuons* theMuons) { 

  if(lpt_mu_n != 0) 
    std::cout << "lpt_mu_n not initialized to 0??" <<std::endl;

  for (int i=0; i<theMuons->mu_n; i++) {
    TLorentzVector * pp4= (TLorentzVector *) theMuons->mu_p4->At(i);
    if(pp4->Pt()<leptons_cut_et_rec) continue;
    if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) continue;
    if(lpt_n==MAX_LEPTONS-1) continue;
    
    new ((*lpt_p4)[lpt_n]) TLorentzVector(*pp4);
    lpt_pdgid[lpt_n] = 13;
    if(theMuons->mu_charge[i]<0) lpt_pdgid[lpt_n]=-13;
    lpt_ind[lpt_n]=i;
    lpt_duplicate[lpt_n]=0;
    
    lpt_n++;
    lpt_mu_n++;
    lpt_emu_n++;
  }
}

void GlobeLeptons::addElectrons(GlobeElectrons * theElectrons) {
  if(lpt_el_n != 0) std::cout << "lpt_el_n not initialized to 0??" <<std::endl;

  Int_t duplicate_n=0;
  Int_t duplicate_electrons[MAX_ELECTRONS];
  for(int i=0; i!=theElectrons->el_n; i++){
    duplicate_electrons[i]=0;
  }
  for (int i=0; i!=theElectrons->el_n; i++) {
    bool found_new_duplicate=false;
    int iscind = theElectrons->el_scind[i];
    for (int ii=i+1; ii!=theElectrons->el_n; ii++) {
      if(theElectrons->el_scind[ii] == iscind) {
        if(duplicate_electrons[i]==0){
          found_new_duplicate=true;
          duplicate_electrons[i]=1+duplicate_n;
        }
        duplicate_electrons[ii]=duplicate_electrons[i];
      }
    }
    if(found_new_duplicate)duplicate_n++;
  }

  // loop for non-duplicate electrons
  for (int i=0; i<theElectrons->el_n; i++) {
    if(duplicate_electrons[i]==0){
      TLorentzVector * pp4= (TLorentzVector *) theElectrons->el_p4->At(i);
      if(pp4->Pt()<leptons_cut_et_rec) continue;
      if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) continue;
      if(lpt_n>=MAX_LEPTONS-1) continue;
      new ((*lpt_p4)[lpt_n]) TLorentzVector(*pp4);
      lpt_pdgid[lpt_n] = 11;
      if(theElectrons->el_charge[i]<0) lpt_pdgid[lpt_n]=-11;
      lpt_ind[lpt_n]=i;
      lpt_duplicate[lpt_n]=0;

      lpt_n++;
      lpt_el_n++;
      lpt_emu_n++;
    }
  }

  // loop for best duplicate electrons (e-over-p closest to 1)
  for(int idup=1;idup!=1+duplicate_n;++idup) {
    TLorentzVector * keeperp4 = new TLorentzVector();;
    Int_t keeperind=-1;
    float besteoverp=1000.;
    for (int i=0; i<theElectrons->el_n; i++) {
      if(duplicate_electrons[i]==idup){
        TLorentzVector * pp4= (TLorentzVector *) theElectrons->el_p4->At(i);
        if(pp4->Pt()<leptons_cut_et_rec) continue;
        if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) continue;
        if(lpt_n>=MAX_LEPTONS-1) continue;

        float thiseoverp = theElectrons->el_eopin[i];
        if(fabs(1.-thiseoverp) < fabs(1.-besteoverp)) {
          besteoverp=thiseoverp;
          keeperp4= pp4;
          keeperind=i;
        }
      }
    }
    if(keeperind == -1)std::cout << "DUPLICATE NOT FOUND!!! IMPOSSIBLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    new ((*lpt_p4)[lpt_n]) TLorentzVector(*keeperp4);
    lpt_pdgid[lpt_n] = 11;
    if(theElectrons->el_charge[keeperind]<0) lpt_pdgid[lpt_n]=-11;
    lpt_ind[lpt_n]=keeperind;
    lpt_duplicate[lpt_n]=idup;
    lpt_n++;
    lpt_el_n++;
    lpt_emu_n++;

  }

  // loop for non-best duplicate electrons
  for(int idup=1;idup!=1+duplicate_n;++idup) {
    //TLorentzVector * keeperp4;
    Int_t keeperind=-1;
    float besteoverp=1000.;
    // first find the best (again)
    for (int i=0; i<theElectrons->el_n; i++) {
      if(duplicate_electrons[i]==idup){
        TLorentzVector * pp4= (TLorentzVector *) theElectrons->el_p4->At(i);
        if(pp4->Pt()<leptons_cut_et_rec) continue;
        if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) continue;
        if(lpt_n>=MAX_LEPTONS-1) continue;

        float thiseoverp = theElectrons->el_eopin[i];
        if(fabs(1.-thiseoverp) < fabs(1.-besteoverp)) {
          besteoverp=thiseoverp;
          //keeperp4= pp4;
          keeperind=i;
        }
      }
    }
    // now fill all other than best
    for (int i=0; i<theElectrons->el_n; i++) {
      if(duplicate_electrons[i]==idup && i != keeperind){
        TLorentzVector * pp4= (TLorentzVector *) theElectrons->el_p4->At(i);
        if(pp4->Pt()<leptons_cut_et_rec) continue;
        if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) continue;
        if(lpt_n>=MAX_LEPTONS-1) continue;

        new ((*lpt_p4)[lpt_n]) TLorentzVector(*pp4);
        lpt_pdgid[lpt_n] = 11;
        if(theElectrons->el_charge[keeperind]<0) lpt_pdgid[lpt_n]=-11;
        lpt_ind[lpt_n]=i;
        lpt_duplicate[lpt_n]=-idup;
        lpt_n++;
        lpt_el_n++;
      }
    }
  }

}

void GlobeLeptons::addPhotons(GlobePhotons * thePhotons) {
  if(lpt_pho_n != 0)
    std::cout << "lpt_pho_n not initialized to 0??" << std::endl;
  // photons

  for (int i=0; i<thePhotons->pho_n; i++) {
    TLorentzVector * pp4= (TLorentzVector *) thePhotons->pho_p4->At(i);
    if(pp4->Pt()<leptons_cut_et_rec) 
      continue;
    if(fabs((float) pp4->Eta())>leptons_cut_eta_rec) 
      continue;
    if(lpt_n>=MAX_LEPTONS-1)
      continue;

    new ((*lpt_p4)[lpt_n]) TLorentzVector(*pp4);
    lpt_pdgid[lpt_n] = 22;
    lpt_ind[lpt_n]=i;
    lpt_duplicate[lpt_n]=0;
    
    lpt_n++;
    lpt_pho_n++;
  } 
}

void GlobeLeptons::fillList(GlobeElectrons * theElectrons, GlobeMuons * theMuons, GlobePhotons * thePhotons) {
  
  lpt_n=0;
  lpt_mu_n=0;
  lpt_el_n=0;
  lpt_emu_n=0;
  lpt_pho_n=0;

  addMuons(theMuons);

  addElectrons(theElectrons);

  addPhotons(thePhotons);

  bool printLeptonCounts=false; 
  if(printLeptonCounts) {
    if(theMuons) std::cout << "mu_n: " << theMuons->mu_n << std::endl;
    if(theElectrons) std::cout << "el_n: " << theElectrons->el_n << std::endl;
    if(thePhotons) std::cout << "pho_n: " << thePhotons->pho_n << std::endl;
    std::cout << "lpt_n: " << lpt_n << std::endl;
    std::cout << "lpt_mu_n: " << lpt_mu_n << std::endl;
    std::cout << "lpt_el_n: " << lpt_el_n << std::endl;
    std::cout << "lpt_emu_n: " << lpt_emu_n << std::endl;
    std::cout << "lpt_pho_n: " << lpt_pho_n << std::endl;
    for(int i=0; i<lpt_n; i++) {
      std::cout << lpt_pdgid[i] << std::endl;
    }
  }
}


