#ifndef GLOBEJETS_H
#define GLOBEJETS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class PFJetIDSelectionFunctor;
class GlobeAnalyzer;

class GlobeJets {
 public:
  
  GlobeJets(const edm::ParameterSet&, const char*);
  virtual ~GlobeJets() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);
 
  // variables
  Int_t jet_n;
  UInt_t jet_nvtx;

  Float_t jet_area[MAX_JETS];
  Float_t jet_emfrac[MAX_JETS];
  Float_t jet_hadfrac[MAX_JETS];

  Float_t jet_pull_dy[MAX_JETS];
  Float_t jet_pull_dphi[MAX_JETS];
  Float_t jet_erescale[MAX_JETS];

  // Jet ID MVA variables
  Float_t jet_dRMean[MAX_JETS];
  Float_t jet_frac01[MAX_JETS];
  Float_t jet_frac02[MAX_JETS];
  Float_t jet_frac03[MAX_JETS];
  Float_t jet_frac04[MAX_JETS];
  Float_t jet_frac05[MAX_JETS];
  Float_t jet_frac06[MAX_JETS];
  Float_t jet_frac07[MAX_JETS];
  Int_t jet_nNeutrals[MAX_JETS];
  Int_t jet_nNeutrals_ptCut[MAX_JETS];
  Float_t jet_beta[MAX_JETS];
  Float_t jet_betaStar[MAX_JETS];
  Float_t jet_dZ[MAX_JETS];
  Int_t jet_nCharged[MAX_JETS];
  Int_t jet_nCharged_QC[MAX_JETS];
  Int_t jet_nCharged_ptCut_QC[MAX_JETS];
  Int_t jet_nCharged_ptCut[MAX_JETS];
  Float_t jet_dR2Mean[MAX_JETS];
  Float_t jet_betaStarClassic[MAX_JETS];
  std::vector<std::vector<float> >* jet_beta_ext;
  std::vector<std::vector<float> >* jet_betaStar_ext;
  std::vector<std::vector<float> >* jet_betaStarClassic_ext;
  
  Bool_t jet_pfloose[MAX_JETS];


  // quark-gluon discrimination:
  Float_t jet_rmsCand[MAX_JETS];
  Float_t jet_ptD[MAX_JETS];
  Float_t jet_axis1[MAX_JETS];
  Float_t jet_axis2[MAX_JETS];
  Float_t jet_pull[MAX_JETS];
  Float_t jet_tana[MAX_JETS];
  Float_t jet_rmsCand_QC[MAX_JETS];
  Float_t jet_ptD_QC[MAX_JETS];
  Float_t jet_axis1_QC[MAX_JETS];
  Float_t jet_axis2_QC[MAX_JETS];
  Float_t jet_pull_QC[MAX_JETS];
  Float_t jet_tana_QC[MAX_JETS];

  Float_t jet_Rchg[MAX_JETS];
  Float_t jet_Rneutral[MAX_JETS];
  Float_t jet_R[MAX_JETS];
  Float_t jet_Rchg_QC[MAX_JETS];

      
  // btags:
  Float_t jet_csvBtag[MAX_JETS];
  Float_t jet_csvMvaBtag[MAX_JETS];
  Float_t jet_jetProbBtag[MAX_JETS];
  Float_t jet_tcheBtag[MAX_JETS];

	// extra info on secondary vertex
	Float_t jet_nSecondaryVertices[MAX_JETS];
	Float_t jet_secVtxPt[MAX_JETS];
	Float_t jet_secVtx3dL[MAX_JETS];
	Float_t jet_secVtx3deL[MAX_JETS];

  Int_t jet_ntk[MAX_JETS];
  Int_t jet_ncalotw[MAX_JETS];
  std::vector<std::vector<unsigned short> >* jet_calotwind;
  std::vector<std::vector<unsigned short> >* jet_tkind;
  std::vector<float * > mvas_;
  std::vector<int * > wp_levels_;
  std::vector<std::vector<std::vector<float> > * > mvas_ext_;
  std::vector<std::vector<std::vector<int> > * > wp_levels_ext_;
  std::vector<PileupJetIdAlgo* > algos_;
    
  TClonesArray *jet_p4;

 private:
  const char* nome;
  GlobeCuts *gCUT;
  edm::InputTag jetColl, calotowerColl, trackColl, jetTkAssColl, vertexColl;
  std::vector<edm::ParameterSet > jetMVAAlgos;
  PFJetIDSelectionFunctor * pfLooseId;
  std::string pfak5corrdata, pfak5corrmc, pfak5corr;
  edm::InputTag bcBColl, bcEColl, tkColl, pfJetColl;
  int debug_level;

  bool doEgammaSummer09Skim;
};

#endif
