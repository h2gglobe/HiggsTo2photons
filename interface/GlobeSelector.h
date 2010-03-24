#ifndef GLOBESELECTOR_H
#define GLOBESELECTOR_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons//interface/Limits.h"

#include <bitset>

class GlobeElectrons;
class GlobeMuons;
class GlobePhotons;
class GlobeGenerator;
class GlobeGenParticles;
class GlobeLeptons;
class GlobeReducedGen;

class GlobeSelector {
 public:
  
  GlobeSelector(const edm::ParameterSet&);
  virtual ~GlobeSelector() {};

  bool selectRecoElectrons(GlobeElectrons*);
  bool selectRecoMuons(GlobeMuons*);
  bool selectRecoPhotons(GlobePhotons*);
  bool selectMCLepPhot(GlobeGenerator*);
  bool selectMCLepPhot(GlobeGenParticles*);
  bool selectGlobal(GlobeLeptons*, GlobeReducedGen*);

#define MAXSELECTORS 5
  
  std::bitset<5> select(GlobeElectrons*, GlobeMuons*, GlobePhotons*, GlobeGenerator*, GlobeLeptons*, GlobeReducedGen*);
std::bitset<5> select(GlobeElectrons*, GlobeMuons*, GlobePhotons*, GlobeGenParticles*, GlobeLeptons*, GlobeReducedGen*);
std::bitset<5> select(GlobeElectrons*, GlobeMuons*, GlobePhotons*);

 private:
  int number_selectors;
  int active_selector;
  int nRecoElectrons[MAXSELECTORS];
  int nRecoMuons[MAXSELECTORS];
  int nRecoPhotons[MAXSELECTORS];
  int nMCLepPhot[MAXSELECTORS];

  double MCetcutel[MAXSELECTORS];
  double MCetcutmu[MAXSELECTORS];
  double MCetcutph[MAXSELECTORS];
  double MCetacutel[MAXSELECTORS];
  double MCetacutmu[MAXSELECTORS];
  double MCetacutph[MAXSELECTORS];

  int nMCel[MAXSELECTORS];
  int nMCmu[MAXSELECTORS];
  int nMCph[MAXSELECTORS];
  int nMCelmu[MAXSELECTORS];
  int nMCelph[MAXSELECTORS];
  int nMCmuph[MAXSELECTORS];
  int nMCelmuph[MAXSELECTORS];
  
  double Recoetcutel[MAXSELECTORS];
  double Recoetcutmu[MAXSELECTORS];
  double Recoetcutph[MAXSELECTORS];
  double Recoetacutel[MAXSELECTORS];
  double Recoetacutmu[MAXSELECTORS];
  double Recoetacutph[MAXSELECTORS];

  int nRecoel[MAXSELECTORS];
  int nRecomu[MAXSELECTORS];
  int nRecoph[MAXSELECTORS];
  int nRecoelmu[MAXSELECTORS];
  int nRecoelph[MAXSELECTORS];
  int nRecomuph[MAXSELECTORS];
  int nRecoelmuph[MAXSELECTORS];


};

#endif
