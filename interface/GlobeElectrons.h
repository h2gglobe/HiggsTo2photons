#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#ifndef GLOBEELECTRONS_H
#define GLOBEELECTRONS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"

//#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
//#include "EGamma/EGammaAnalysisTools/interface/ElectronEnergyRegressionEvaluate.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class ElectronMVAEstimator;

class GlobeElectrons {
 public:
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

  GlobeElectrons(const edm::ParameterSet&, const char* n = "std");
  virtual ~GlobeElectrons();

  void defineBranch(TTree* tree);
  bool analyze(const edm::Event&, const edm::EventSetup&);
  bool analyze_pf(const edm::Event&, const edm::EventSetup&);
  void initialize_branches(int electron_number);
  float hoeCalculator(const reco::BasicCluster*, const CaloGeometry&,
                      const edm::Event&, const edm::EventSetup&);
  std::map<DetId, EcalRecHit> rechits_map_;

  EGammaMvaEleEstimator* myMVANonTrig;
  EGammaMvaEleEstimator* myMVATrig;
  std::vector<std::string> myManualCatWeightsNonTrig;
  std::vector<std::string> myManualCatWeightsTrig;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

  //bool identify(const reco::GsfElectronRef electron, int type);
  //bool st_identify(const reco::GsfElectronRef electron, int type);
  //int classify(const reco::GsfElectronRef electron);

  std::pair<unsigned int, float> sharedHits(const reco::Track& trackA, const reco::Track& trackB);

  // variables
  Int_t el_n;
  Float_t el_enearbcopin[MAX_ELECTRONS];
  Float_t el_eseedopout[MAX_ELECTRONS];
  Float_t el_eopin[MAX_ELECTRONS];
  Float_t el_pout[MAX_ELECTRONS];
  Float_t el_pin[MAX_ELECTRONS];
  Float_t el_fbrem[MAX_ELECTRONS];
  Int_t el_nbrem[MAX_ELECTRONS];
  Int_t el_1pxb[MAX_ELECTRONS];
  Int_t el_1pxf[MAX_ELECTRONS];

  Float_t el_hoe[MAX_ELECTRONS];
  Float_t el_hoed1[MAX_ELECTRONS];
  Float_t el_hoed2[MAX_ELECTRONS];

  Float_t el_hoebc[MAX_ELECTRONS];
  Float_t el_hoebcd1[MAX_ELECTRONS];
  Float_t el_hoebcd2[MAX_ELECTRONS];

  Float_t el_detain[MAX_ELECTRONS];
  Float_t el_dphiin[MAX_ELECTRONS];
  Float_t el_detaout[MAX_ELECTRONS];
  Float_t el_dphiout[MAX_ELECTRONS];
  Float_t el_z0[MAX_ELECTRONS];
  Float_t el_d0[MAX_ELECTRONS];
  Float_t el_z0err[MAX_ELECTRONS];
  Float_t el_d0err[MAX_ELECTRONS];
  Float_t el_chi2[MAX_ELECTRONS];
  Float_t el_dof[MAX_ELECTRONS];
  Float_t el_e1x5[MAX_ELECTRONS];
  Float_t el_e5x5[MAX_ELECTRONS];
  Float_t el_sipip[MAX_ELECTRONS];
  Float_t el_sieie[MAX_ELECTRONS];
  Float_t el_sieip[MAX_ELECTRONS];
  Float_t el_sieiesc[MAX_ELECTRONS];
  Float_t el_eseffsixix[MAX_ELECTRONS];
  Float_t el_eseffsiyiy[MAX_ELECTRONS];
  Float_t el_e2x5[MAX_ELECTRONS];
  Float_t el_esc[MAX_ELECTRONS];
  Float_t el_eseed[MAX_ELECTRONS];
  Float_t el_eseedopin[MAX_ELECTRONS];
  Float_t el_qoverperr[MAX_ELECTRONS];
  Float_t el_pterr[MAX_ELECTRONS];
  Float_t el_etaerr[MAX_ELECTRONS];
  Float_t el_phierr[MAX_ELECTRONS];
  Int_t el_class[MAX_ELECTRONS];
  Int_t el_charge[MAX_ELECTRONS];
  Int_t el_ch_gsf[MAX_ELECTRONS];
  Int_t el_ch_scpix[MAX_ELECTRONS];
  Int_t el_losthits[MAX_ELECTRONS];
  Int_t el_validhits[MAX_ELECTRONS];
  Int_t el_hp_expin[MAX_ELECTRONS];
  Int_t el_hp_expout[MAX_ELECTRONS];
  Float_t el_must[MAX_ELECTRONS];
  Int_t el_mustnc[MAX_ELECTRONS];
  Float_t el_1oe_1op[MAX_ELECTRONS];

  Int_t el_scind[MAX_ELECTRONS];
  Int_t el_crack[MAX_ELECTRONS];
  Int_t el_tkind[MAX_ELECTRONS];
  Int_t el_conv[MAX_ELECTRONS];
  
  Int_t el_nambtk[MAX_ELECTRONS];

  Float_t el_pfiso_myneutral03[MAX_ELECTRONS];
  Float_t el_pfiso_mycharged03[MAX_ELECTRONS];
  Float_t el_pfiso_myphoton03[MAX_ELECTRONS];
  Float_t el_pfiso_myneutral04[MAX_ELECTRONS];
  Float_t el_pfiso_mycharged04[MAX_ELECTRONS];
  Float_t el_pfiso_myphoton04[MAX_ELECTRONS];
  Float_t el_pfiso_neutral[MAX_ELECTRONS];
  Float_t el_pfiso_charged[MAX_ELECTRONS];
  Float_t el_pfiso_photon[MAX_ELECTRONS];
  Float_t el_tkiso03[MAX_ELECTRONS];
  Float_t el_ecaliso03[MAX_ELECTRONS];
  Float_t el_hcaliso03[MAX_ELECTRONS];
  Float_t el_hcalsolidiso03[MAX_ELECTRONS];
  Float_t el_tkiso04[MAX_ELECTRONS];
  Float_t el_ecaliso04[MAX_ELECTRONS];
  Float_t el_hcaliso04[MAX_ELECTRONS];
  Float_t el_hcalsolidiso04[MAX_ELECTRONS];
  Float_t el_hcalbciso03[MAX_ELECTRONS];
  Float_t el_hcalbciso04[MAX_ELECTRONS];

  Float_t el_r9[MAX_ELECTRONS];
  Float_t el_gsfchi2[MAX_ELECTRONS];
  Float_t el_ip3d_err[MAX_ELECTRONS];
  Float_t el_ip3d[MAX_ELECTRONS];
  Float_t el_ip3d_sig[MAX_ELECTRONS];
  Float_t el_sc_time[MAX_ELECTRONS];

  //Float_t el_mva[MAX_ELECTRONS];  
  //Float_t el_mva_noiso[MAX_ELECTRONS];
  Float_t el_mva_nontrig[MAX_ELECTRONS];  
  Float_t el_mva_trig[MAX_ELECTRONS];  
  Bool_t el_ecaldrv[MAX_ELECTRONS];
  Bool_t el_tkdrv[MAX_ELECTRONS];
  Float_t el_ip_ctf[MAX_ELECTRONS];
  Float_t el_ip_gsf[MAX_ELECTRONS];
  Float_t el_dist[MAX_ELECTRONS];
  Float_t el_dcot[MAX_ELECTRONS];
  
  Float_t el_regr_energy[MAX_ELECTRONS];
  Float_t el_regr_energyerr[MAX_ELECTRONS];
  Float_t el_corr_energy[MAX_ELECTRONS];
  Float_t el_corr_energyerr[MAX_ELECTRONS];
  Float_t el_calib_energy[MAX_ELECTRONS];
  Float_t el_calib_energyerr[MAX_ELECTRONS];
  Float_t el_eleopout[MAX_ELECTRONS];
  Float_t el_detaeleout[MAX_ELECTRONS];
  Int_t el_kfhits[MAX_ELECTRONS];
  Float_t  el_kfchi2[MAX_ELECTRONS];

  Float_t el_psenergy[MAX_ELECTRONS];  
  Int_t el_passmvapresel[MAX_ELECTRONS];
  Int_t el_passcutpresel[MAX_ELECTRONS];
  Float_t el_psenergypf[MAX_ELECTRONS];  
  Int_t el_nbrempf[MAX_ELECTRONS];
  Float_t el_eseedpf[MAX_ELECTRONS];
  Float_t el_epf[MAX_ELECTRONS];
  Float_t el_psly1[MAX_ELECTRONS];
  Float_t el_psly2[MAX_ELECTRONS];
  Int_t el_psnstriply1[MAX_ELECTRONS];
  Int_t el_psnstriply2[MAX_ELECTRONS];
  
  Float_t el_D0Vtx[MAX_ELECTRONS][100];
  Float_t el_DZVtx[MAX_ELECTRONS][100];
  Float_t el_conv_vtxProb[MAX_ELECTRONS];

  /** corresponds to the el_XXX_catbased variable in the output tree.

      - first index is the electron index
      - second index corresponds to the eID label (type of electron id)
        as specified in the eIDLabels parameter of GlobeAnalyzer
        (e.g. eidLoose, eidTight)
      - the value is a bit pattern, typically 

           0 - no cut passed
           1 - eID cuts passed
           2 - iso cuts passed
           4 - conversion rejection
           8 - ip cut 
       
      See also https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCategoryBasedElectronID#How_to_use_it_in_CMSSW
  */ 
  std::vector<std::vector<int> >* el_catbased;

  TClonesArray *el_sc;
  TClonesArray *el_p4;
  TClonesArray *el_p4_corr;
  TClonesArray *el_momvtx;
  TClonesArray *el_momvtxconst;
  TClonesArray *el_momcalo;
  TClonesArray *el_momout;
  TClonesArray *el_posvtx;
  TClonesArray *el_poscalo;

 private:
  const char* nome; 
  bool doAodSim;
  GlobeCuts *gCUT;
  edm::InputTag electronColl, trackColl, trackColl2, vertexColl, beamSpotColl, conversionColl, pfColl;
  GlobeEcalClusters *gES;
  std::vector<edm::InputTag> eIDLabels;
  int debug_level;

  // SUPER CLUSTERS
  edm::InputTag hybridSuperClusterColl;
  edm::InputTag endcapSuperClusterColl;
  edm::InputTag rhoCollection;  
  edm::InputTag caloTowerColl;
  edm::InputTag ecalHitEBColl;
  edm::InputTag ecalHitEEColl;
  edm::InputTag ecalHitESColl;
  edm::InputTag hcalHitColl;
  //edm::FileInPath mvaWeightFile;  
  std::vector<std::string> mvaNonTrigWeightFiles;  
  std::vector<std::string> mvaTrigWeightFiles;  
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;

  EGEnergyCorrector ecorr_;
  bool energyCorrectionsFromDB;
  std::string energyRegFilename;
  std::string regressionVersion;
  //std::string eleRegressionFilename;
  //Int_t eleRegressionType;

  //ElectronEnergyRegressionEvaluate* eleRegression;

  const TransientTrackBuilder* transientTrackBuilder;
};

#endif
