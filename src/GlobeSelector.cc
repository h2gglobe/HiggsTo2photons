#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSelector.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenerator.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeReducedGen.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"

GlobeSelector::GlobeSelector(const edm::ParameterSet& iConfig) {

  number_selectors = iConfig.getParameter<int>("SelectorNumber");

  if (number_selectors > MAXSELECTORS) {
    std::cout << "TOO MANY SELECTORS, ONLY 5 ALLOWED !!! " << std::endl;
    number_selectors = MAXSELECTORS;
  }

  //set them to a large number (I don't know if it is needed (Marco))
  for(int i = 0; i<MAXSELECTORS; ++i) {
    nRecoElectrons[i]=1000;
    nRecoMuons[i]=1000;
    nRecoPhotons[i]=1000;
    nMCLepPhot[i]=1000;

    nMCel[i]=1000;
    nMCmu[i]=1000;
    nMCph[i]=1000;
    nMCelmu[i]=1000;
    nMCelph[i]=1000;
    nMCmuph[i]=1000;
    nMCelmuph[i]=1000;

    nRecoel[i]=1000;
    nRecomu[i]=1000;
    nRecoph[i]=1000;
    nRecoelmu[i]=1000;
    nRecoelph[i]=1000;
    nRecomuph[i]=1000;
    nRecoelmuph[i]=1000;

  }

  for(int i = 0; i<number_selectors; ++i) {
    char a[100];
    sprintf(a, "SelectorCuts%d", i);
    edm::ParameterSet psetSelector = iConfig.getParameter<edm::ParameterSet>(a);
    nRecoElectrons[i] = psetSelector.getParameter<int>("RecoElectronNumber"); 
    nRecoMuons[i] = psetSelector.getParameter<int>("RecoMuonNumber"); 
    nRecoPhotons[i] = psetSelector.getParameter<int>("RecoPhotonNumber"); 
    nMCLepPhot[i] = psetSelector.getParameter<int>("MCLepPhotNumber"); 

    MCetcutel[i] = psetSelector.getParameter<double>("MCetcutel"); 
    MCetcutmu[i] = psetSelector.getParameter<double>("MCetcutmu"); 
    MCetcutph[i] = psetSelector.getParameter<double>("MCetcutph"); 
    MCetacutel[i] = psetSelector.getParameter<double>("MCetacutel"); 
    MCetacutmu[i] = psetSelector.getParameter<double>("MCetacutmu"); 
    MCetacutph[i] = psetSelector.getParameter<double>("MCetacutph"); 

    nMCel[i] = psetSelector.getParameter<int>("nMCel"); 
    nMCmu[i] = psetSelector.getParameter<int>("nMCmu"); 
    nMCph[i] = psetSelector.getParameter<int>("nMCph"); 
    nMCelmu[i] = psetSelector.getParameter<int>("nMCelmu"); 
    nMCelph[i] = psetSelector.getParameter<int>("nMCelph"); 
    nMCmuph[i] = psetSelector.getParameter<int>("nMCmuph"); 
    nMCelmuph[i] = psetSelector.getParameter<int>("nMCelmuph"); 

    Recoetcutel[i] = psetSelector.getParameter<double>("Recoetcutel"); 
    Recoetcutmu[i] = psetSelector.getParameter<double>("Recoetcutmu"); 
    Recoetcutph[i] = psetSelector.getParameter<double>("Recoetcutph"); 
    Recoetacutel[i] = psetSelector.getParameter<double>("Recoetacutel"); 
    Recoetacutmu[i] = psetSelector.getParameter<double>("Recoetacutmu"); 
    Recoetacutph[i] = psetSelector.getParameter<double>("Recoetacutph"); 

    nRecoel[i] = psetSelector.getParameter<int>("nRecoel"); 
    nRecomu[i] = psetSelector.getParameter<int>("nRecomu"); 
    nRecoph[i] = psetSelector.getParameter<int>("nRecoph"); 
    nRecoelmu[i] = psetSelector.getParameter<int>("nRecoelmu"); 
    nRecoelph[i] = psetSelector.getParameter<int>("nRecoelph"); 
    nRecomuph[i] = psetSelector.getParameter<int>("nRecomuph"); 
    nRecoelmuph[i] = psetSelector.getParameter<int>("nRecoelmuph"); 

  }
}

bool GlobeSelector::selectGlobal(GlobeLeptons* lept, GlobeReducedGen* redgen) {
  //std::cout<<" inside GlobeSelector::selectGlobal "<<std::endl;
  int nrecopassel=0;
  int nrecopassmu=0;
  int nrecopassph=0;
  int nrecopasselmu=0;
  int nrecopasselph=0;
  int nrecopassmuph=0;
  int nrecopasselmuph=0;

  int nmcpassel=0;
  int nmcpassmu=0;
  int nmcpassph=0;
  int nmcpasselmu=0;
  int nmcpasselph=0;
  int nmcpassmuph=0;
  int nmcpasselmuph=0;

  //missing for now duplicates between electrons and photons
  if(lept) { 
    for (int i = 0 ; i < lept->lpt_n; i++){ 
      TLorentzVector * p4 = (TLorentzVector *) lept->lpt_p4->At(i);
      //std::cout<<"REC LEPTON "<<i<<" "<<lept->lpt_pdgid[i]<<" dupl "<<lept->lpt_duplicate[i]<<" "<<p4->Eta()<<" "<<p4->Et()<<std::endl;

      //electrons
      if(abs(lept->lpt_pdgid[i])==11) {   
        if(lept->lpt_duplicate[i]>=0) {
          if( fabs(p4->Eta()) < Recoetacutel[active_selector] && p4->Et() > Recoetcutel[active_selector]) {
            nrecopassel++;
            nrecopasselmu++;
            nrecopasselph++;
            nrecopasselmuph++;
          }
        }
      }
      //muons
      if(abs(lept->lpt_pdgid[i])==13) {   
        if(lept->lpt_duplicate[i]>=0) {
          if( fabs(p4->Eta()) < Recoetacutmu[active_selector] && p4->Et() > Recoetcutmu[active_selector]) {
            nrecopassmu++;
            nrecopasselmu++;
            nrecopassmuph++;
            nrecopasselmuph++;
          }
        }
      }
      //photons
      if(abs(lept->lpt_pdgid[i])==22) {   
        if(lept->lpt_duplicate[i]>=0) {
          if( fabs(p4->Eta()) < Recoetacutph[active_selector] && p4->Et() > Recoetcutph[active_selector]) {
            nrecopassph++;
            nrecopasselph++;
            nrecopassmuph++;
            nrecopasselmuph++;
          }
        }
      }
    }
  }
  if(redgen) { 
    for (int i = 0 ; i < redgen->lptgen_n; i++){ 
      TLorentzVector * p4 = (TLorentzVector *) redgen->lptgen_p4->At(i);
      //std::cout<<"GEN LEPTON "<<i<<" "<<redgen->lptgen_pdgid[i]<<" "<<p4->Eta()<<" "<<p4->Et()<<std::endl;

      //electrons
      if(abs(redgen->lptgen_pdgid[i])==11) {   
        if( fabs(p4->Eta()) < Recoetacutel[active_selector] && p4->Et() > Recoetcutel[active_selector]) {
          nmcpassel++;
          nmcpasselmu++;
          nmcpasselph++;
          nmcpasselmuph++;
        }
      }
      //muons
      if(abs(redgen->lptgen_pdgid[i])==13) {   
        if( fabs(p4->Eta()) < Recoetacutmu[active_selector] && p4->Et() > Recoetcutmu[active_selector]) {
          nmcpassmu++;
          nmcpasselmu++;
          nmcpassmuph++;
          nmcpasselmuph++;
        }
      }
      //photons
      if(abs(redgen->lptgen_pdgid[i])==22) {   
        if( fabs(p4->Eta()) < Recoetacutph[active_selector] && p4->Et() > Recoetcutph[active_selector]) {
          nmcpassph++;
          nmcpasselph++;
          nmcpassmuph++;
          nmcpasselmuph++;
        }
      }
    }
  }


  //std::cout<<active_selector<<" "<<nmcpassel<<nmcpassmu<<nmcpassph<<nmcpasselmu<<nmcpasselph<<nmcpassmuph<<nmcpasselmuph<<std::endl;
  //std::cout<<" nMCel "<<nMCel[active_selector]<<std::endl;

  if(nmcpassel>=nMCel[active_selector]) return 1;
  if(nmcpassmu>=nMCmu[active_selector]) return 1;
  if(nmcpassph>=nMCph[active_selector]) return 1;
  if(nmcpasselmu>=nMCelmu[active_selector]) return 1;
  if(nmcpasselph>=nMCelph[active_selector]) return 1;
  if(nmcpassmuph>=nMCmuph[active_selector]) return 1;
  if(nmcpasselmuph>=nMCelmuph[active_selector]) return 1;

  if(nrecopassel>=nRecoel[active_selector]) return 1;
  if(nrecopassmu>=nRecomu[active_selector]) return 1;
  if(nrecopassph>=nRecoph[active_selector]) return 1;
  if(nrecopasselmu>=nRecoelmu[active_selector]) return 1;
  if(nrecopasselph>=nRecoelph[active_selector]) return 1;
  if(nrecopassmuph>=nRecomuph[active_selector]) return 1;
  if(nrecopasselmuph>=nRecoelmuph[active_selector]) return 1;

  return 0;
}

bool GlobeSelector::selectRecoElectrons(GlobeElectrons* el) {
  int nel=0;
  if(el) for (int i = 0 ; i < el->el_n; i++){
    TLorentzVector * elp4 = (TLorentzVector *) el->el_p4->At(i);
    // std::cout<<"El energy = "<<elp4->E()<<"\t id :"<<el->el_robust[i]<<std::endl;
    if( fabs(elp4->Eta()) < 3. && elp4->Et() > 5.)// && el->el_roloose[i] == 1)
    {
      nel++;
    }
  }
  if (nel>= nRecoElectrons[active_selector]) return 1;

  return 0;
  //  return (el->el_n >= nRecoElectrons[active_selector]); 
}

bool GlobeSelector::selectRecoMuons(GlobeMuons* mu) {
  int nmu=0;
  if(mu) for (int i = 0 ; i < mu->mu_n; i++){
    TLorentzVector * mup4 = (TLorentzVector *) mu->mu_p4->At(i);
    if( fabs(mup4->Eta()) < 3. && mup4->Et() > 5. )
    {
      nmu++;
    }
  }
  if (nmu >= nRecoMuons[active_selector]) return 1;

  return 0;
  //  return (mu->mu_n >= nRecoMuons[active_selector]); 
}

bool GlobeSelector::selectRecoPhotons(GlobePhotons* phot){
  int nph=0;
  if(phot) for (int i = 0 ; i < phot->pho_n; i++){
    TLorentzVector * photp4 = (TLorentzVector *) phot->pho_p4->At(i);
    if( fabs(photp4->Eta()) < 3. && photp4->Et() > 10. )
    {
      nph++;
    }
  }
  if (nph >= nRecoPhotons[active_selector]) return 1;
  return 0;
  //  return (phot->pho_n >= nRecoPhotons[active_selector]); 
}

bool GlobeSelector::selectMCLepPhot(GlobeGenerator* gen){
  int nlptph=0;
  if(gen) for (int i = 0 ; i < gen->gen_n; i++) if(gen->gen_status[i] == 1) {
    TLorentzVector * genp4 = (TLorentzVector *) gen->gen_p4->At(i);
    if( (fabs(genp4->Eta()) < 3. && genp4->Et() > 5.) && (fabs(gen->gen_pdgid[i]) == 11 || fabs(gen->gen_pdgid[i]) == 13 ))
    {
      nlptph++;
      //return 1;
    }
    if( fabs(genp4->Eta()) < 3. && genp4->Et() > 10. && fabs(gen->gen_pdgid[i]) == 22 )
    {
      nlptph++;
      //return 1;
    }
  }
  if (nlptph >= nMCLepPhot[active_selector]) return 1;
  return 0;
  //return (gen->gen_n >= nMCLepPhot[active_selector]); 
}


bool GlobeSelector::selectMCLepPhot(GlobeGenParticles* gp){
  int nlptph=0;
  if(gp) for (int i = 0 ; i < gp->gp_n; i++) if(gp->gp_status[i] == 1) {
    TLorentzVector * gpp4 = (TLorentzVector *) gp->gp_p4->At(i);
    if( (fabs(gpp4->Eta()) < 3. && gpp4->Et() > 5.) && (fabs(gp->gp_pdgid[i]) == 11 || fabs(gp->gp_pdgid[i]) == 13 ))
    {
      nlptph++;
      //return 1;
    }
    if( fabs(gpp4->Eta()) < 3. && gpp4->Et() > 10. && fabs(gp->gp_pdgid[i]) == 22 )
    {
      nlptph++;
      //return 1;
    }
  }
  if (nlptph >= nMCLepPhot[active_selector]) return 1;
  return 0;
  //return (gp->gp_n >= nMCLepPhot[active_selector]); 
}


std::bitset<5> GlobeSelector::select(GlobeElectrons* el, GlobeMuons* mu, GlobePhotons* phot, GlobeGenerator* gen, GlobeLeptons* lept, GlobeReducedGen* redgen) {
  //std::cout<<" inside GlobeSelector::select number_selectors "<<number_selectors<<std::endl;

  if (number_selectors == 0)
    return true;

  std::bitset<5> bits;
  for(active_selector=0; active_selector<number_selectors; ++active_selector) {
    if (selectRecoElectrons(el) || selectRecoMuons(mu) || selectRecoPhotons(phot) || selectMCLepPhot(gen) || selectGlobal(lept, redgen) ) {
      bits[active_selector] = 1;
      //std::cout<<"*** Event selected by Selection "<<active_selector<<std::endl;
    }
    else {
      //std::cout<<"<<< Event did not pass Selection "<<active_selector<<std::endl;
    }
  }
  return bits;
}


std::bitset<5> GlobeSelector::select(GlobeElectrons* el, GlobeMuons* mu, GlobePhotons* phot, GlobeGenParticles* gp, GlobeLeptons* lept, GlobeReducedGen* redgen) {
  //std::cout<<" inside GlobeSelector::select number_selectors "<<number_selectors<<std::endl;

  if (number_selectors == 0)
    return true;

  std::bitset<5> bits;
  for(active_selector=0; active_selector<number_selectors; ++active_selector) {
    if (selectRecoElectrons(el) || selectRecoMuons(mu) || selectRecoPhotons(phot) || selectMCLepPhot(gp) || selectGlobal(lept, redgen) ) {
      bits[active_selector] = 1;
      //std::cout<<"*** Event selected by Selection "<<active_selector<<std::endl;
    }
    else {
      //std::cout<<"<<< Event did not pass Selection "<<active_selector<<std::endl;
    }
  }
  return bits;
}


std::bitset<5> GlobeSelector::select(GlobeElectrons* el, GlobeMuons* mu, GlobePhotons* phot) {

  if (number_selectors == 0)
    return true;

  std::bitset<5> bits;
  for(active_selector=0; active_selector<number_selectors; ++active_selector) {
    if (selectRecoElectrons(el) || selectRecoMuons(mu) || selectRecoPhotons(phot)) {
      bits[active_selector] = 1;
      //std::cout<<"*** Event selected by Selection "<<active_selector<<std::endl;
    }
    else {
      //std::cout<<"<<< Event did not pass Selection "<<active_selector<<std::endl;
    }
  }
  return bits;
}




