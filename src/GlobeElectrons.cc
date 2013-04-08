#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/PFIsolation.h"

#include <iostream>

GlobeElectrons::GlobeElectrons(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  char a[100];
  sprintf(a, "ElectronColl_%s", nome);
  electronColl = iConfig.getParameter<edm::InputTag>(a);
  debug_level = iConfig.getParameter<int>("Debug_Level");
  doAodSim = iConfig.getParameter<bool>("doAodSim");

  conversionColl = iConfig.getParameter<edm::InputTag>("ConvertedPhotonColl");
  beamSpotColl = iConfig.getParameter<edm::InputTag>("BeamSpot");
  trackColl = iConfig.getParameter<edm::InputTag>("TrackColl");
  trackColl2 = iConfig.getParameter<edm::InputTag>("TrackColl3");
  vertexColl = iConfig.getParameter<edm::InputTag>("VertexColl_std");
  pfColl = iConfig.getParameter<edm::InputTag>("PFCandidateColl");

  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");
  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");
  ecalHitESColl = iConfig.getParameter<edm::InputTag>("EcalHitESColl");
  hcalHitColl = iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");

  eIDLabels = iConfig.getParameter<std::vector<edm::InputTag> >("eIDLabels");

  applyEnergyCorrection = iConfig.getParameter<bool>("applyEnergyCorrection");

  mvaNonTrigWeightFiles = iConfig.getParameter<std::vector<std::string> >("electronNonTrigMVAWeightFileNames");
  mvaTrigWeightFiles    = iConfig.getParameter<std::vector<std::string> >("electronTrigMVAWeightFileNames");
  
  for(unsigned int j=0; j<mvaTrigWeightFiles.size(); j++)
    myManualCatWeightsTrig.push_back(edm::FileInPath(mvaTrigWeightFiles[j]).fullPath());

  for(unsigned int j=0; j<mvaNonTrigWeightFiles.size(); j++)
    myManualCatWeightsNonTrig.push_back(edm::FileInPath(mvaNonTrigWeightFiles[j]).fullPath());

  myMVANonTrig = new EGammaMvaEleEstimator();
  myMVANonTrig->initialize("BDT",
           EGammaMvaEleEstimator::kNonTrig,
           true, // use manual cat
           myManualCatWeightsNonTrig);

  myMVATrig = new EGammaMvaEleEstimator();
  myMVATrig->initialize("BDT",
           EGammaMvaEleEstimator::kTrig,
           true, // use manual cat
           myManualCatWeightsTrig);

  inputTagIsoValElectronsPFId_   = iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");   

  energyCorrectionsFromDB = iConfig.getParameter<bool> ("energyCorrectionsFromDB"); 
  energyRegFilename       = iConfig.getParameter<std::string> ("energyCorrectionsFileNameEle");  
  regressionVersion       = iConfig.getParameter<std::string> ("energyCorrectionsVersion");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
  gES  = new GlobeEcalClusters(iConfig);
}

void GlobeElectrons::defineBranch(TTree* tree) {

  el_sc = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_p4 = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_momvtx = new TClonesArray("TVector3", MAX_ELECTRONS);  
  el_momvtxconst = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momcalo = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momout = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_posvtx = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_poscalo = new TClonesArray("TVector3", MAX_ELECTRONS);

  el_catbased = new std::vector<std::vector<int> >;
  el_schits = new std::vector<std::vector<UInt_t> >;
  el_bchits = new std::vector<std::vector<UInt_t> >;

  char a1[50], a2[50];
  
  sprintf(a1, "el_%s_n", nome);
  sprintf(a2, "el_%s_n/I", nome);
  tree->Branch(a1, &el_n, a2);
  
  sprintf(a1, "el_%s_sc", nome);
  tree->Branch(a1, "TClonesArray", &el_sc, 32000, 0);
  
  sprintf(a1, "el_%s_p4", nome);
  tree->Branch(a1, "TClonesArray", &el_p4, 32000, 0);

  sprintf(a1, "el_%s_momvtx", nome);
  tree->Branch(a1, "TClonesArray", &el_momvtx, 32000, 0);

  sprintf(a1, "el_%s_momvtxconst", nome);
  tree->Branch(a1, "TClonesArray", &el_momvtxconst, 32000, 0);

  sprintf(a1, "el_%s_momcalo", nome);
  tree->Branch(a1, "TClonesArray", &el_momcalo, 32000, 0);

  sprintf(a1, "el_%s_momout", nome);
  tree->Branch(a1, "TClonesArray", &el_momout, 32000, 0);
  
  sprintf(a1, "el_%s_posvtx", nome);
  tree->Branch(a1, "TClonesArray", &el_posvtx, 32000, 0);

  sprintf(a1, "el_%s_poscalo", nome);
  tree->Branch(a1, "TClonesArray", &el_poscalo, 32000, 0);
  
  sprintf(a1, "el_%s_eopin", nome);
  sprintf(a2, "el_%s_eopin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eopin, a2);
  
  sprintf(a1, "el_%s_eseedopout", nome);
  sprintf(a2, "el_%s_eseedopout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseedopout, a2);
  
  sprintf(a1, "el_%s_pout", nome);
  sprintf(a2, "el_%s_pout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pout, a2);
  
  sprintf(a1, "el_%s_pin", nome);
  sprintf(a2, "el_%s_pin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pin, a2);
  
  sprintf(a1, "el_%s_e1x5", nome);
  sprintf(a2, "el_%s_e1x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e1x5, a2);
  
  sprintf(a1, "el_%s_e5x5", nome);
  sprintf(a2, "el_%s_e5x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e5x5, a2);

  sprintf(a1, "el_%s_e2x5", nome);
  sprintf(a2, "el_%s_e2x5[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_e2x5, a2);
 
  sprintf(a1, "el_%s_sipip", nome);
  sprintf(a2, "el_%s_sipip[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_sipip, a2);
  
  sprintf(a1, "el_%s_sieie", nome);
  sprintf(a2, "el_%s_sieie[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_sieie, a2);

  sprintf(a1, "el_%s_sieiesc", nome);
  sprintf(a2, "el_%s_sieiesc[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_sieiesc, a2);

  sprintf(a1, "el_%s_eseffsixix", nome);
  sprintf(a2, "el_%s_eseffsixix[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseffsixix, a2);

  sprintf(a1, "el_%s_eseffsiyiy", nome);
  sprintf(a2, "el_%s_eseffsiyiy[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseffsiyiy, a2);

  sprintf(a1, "el_%s_eseedopin", nome);
  sprintf(a2, "el_%s_eseedopin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseedopin, a2);
  
  sprintf(a1, "el_%s_fbrem", nome);
  sprintf(a2, "el_%s_fbrem[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_fbrem, a2);

  sprintf(a1, "el_%s_nbrem", nome);
  sprintf(a2, "el_%s_nbrem[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_nbrem, a2);

  sprintf(a1, "el_%s_hoe", nome);
  sprintf(a2, "el_%s_hoe[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoe, a2);

  sprintf(a1, "el_%s_hoed1", nome);
  sprintf(a2, "el_%s_hoed1[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoed1, a2);

  sprintf(a1, "el_%s_hoed2", nome);
  sprintf(a2, "el_%s_hoed2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hoed2, a2);

  //sprintf(a1, "el_%s_hoebc", nome);
  //sprintf(a2, "el_%s_hoebc[el_%s_n]/F", nome, nome);
  //tree->Branch(a1, &el_hoebc, a2);

  //sprintf(a1, "el_%s_hoebcd1", nome);
  //sprintf(a2, "el_%s_hoebcd1[el_%s_n]/F", nome, nome);
  //tree->Branch(a1, &el_hoebcd1, a2);

  //sprintf(a1, "el_%s_hoebcd2", nome);
  //sprintf(a2, "el_%s_hoebcd2[el_%s_n]/F", nome, nome);
  //tree->Branch(a1, &el_hoebcd2, a2);
  
  sprintf(a1, "el_%s_detain", nome);
  sprintf(a2, "el_%s_detain[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_detain, a2);
  
  sprintf(a1, "el_%s_dphiin", nome);
  sprintf(a2, "el_%s_dphiin[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dphiin, a2);
  
  sprintf(a1, "el_%s_detaout", nome);
  sprintf(a2, "el_%s_detaout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_detaout, a2);
  
  sprintf(a1, "el_%s_dphiout", nome);
  sprintf(a2, "el_%s_dphiout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dphiout, a2);
   
  sprintf(a1, "el_%s_class", nome);
  sprintf(a2, "el_%s_class[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_class, a2);
 
  sprintf(a1, "el_%s_crack", nome);
  sprintf(a2, "el_%s_crack[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_crack, a2);
   
  sprintf(a1, "el_%s_nambtk", nome);
  sprintf(a2, "el_%s_nambtk[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_nambtk, a2);

  sprintf(a1, "el_%s_scind", nome);
  sprintf(a2, "el_%s_scind[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_scind, a2);
  
  sprintf(a1, "el_%s_z0", nome);
  sprintf(a2, "el_%s_z0[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_z0, a2);
  
  sprintf(a1, "el_%s_d0", nome);
  sprintf(a2, "el_%s_d0[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_d0, a2);

  sprintf(a1, "el_%s_chi2", nome);
  sprintf(a2, "el_%s_chi2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_chi2, a2);
  
  sprintf(a1, "el_%s_mva", nome);
  sprintf(a2, "el_%s_mva[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_mva, a2);
  
  sprintf(a1, "el_%s_ch_gsf", nome);
  sprintf(a2, "el_%s_ch_gsf[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ch_gsf, a2);

  sprintf(a1, "el_%s_ch_scpix", nome);
  sprintf(a2, "el_%s_ch_scpix[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ch_scpix, a2);

  sprintf(a1, "el_%s_charge", nome);
  sprintf(a2, "el_%s_charge[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_charge, a2);
  
  sprintf(a1, "el_%s_losthits", nome);
  sprintf(a2, "el_%s_losthits[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_losthits, a2);
  
  sprintf(a1, "el_%s_validhits", nome);
  sprintf(a2, "el_%s_validhits[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_validhits, a2);

  sprintf(a1, "el_%s_hp_expin", nome);
  sprintf(a2, "el_%s_hp_expin[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_hp_expin, a2);

  sprintf(a1, "el_%s_hp_expout", nome);
  sprintf(a2, "el_%s_hp_expout[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_hp_expout, a2);

  sprintf(a1, "el_%s_catbased", nome);
  tree->Branch(a1, "std::vector<std::vector<int> >", &el_catbased);

  sprintf(a1, "el_%s_schits", nome);
  tree->Branch(a1, "std::vector<std::vector<UInt_t> >", &el_schits);

  sprintf(a1, "el_%s_bchits", nome);
  tree->Branch(a1, "std::vector<std::vector<UInt_t> >", &el_bchits);

  sprintf(a1, "el_%s_tkind", nome);
  sprintf(a2, "el_%s_tkind[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_tkind, a2); 

  sprintf(a1, "el_%s_pfiso_myneutral03", nome);
  sprintf(a2, "el_%s_pfiso_myneutral03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_myneutral03, a2);

  sprintf(a1, "el_%s_pfiso_mycharged03", nome);
  sprintf(a2, "el_%s_pfiso_mycharged03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_mycharged03, a2);
  
  sprintf(a1, "el_%s_pfiso_myphoton03", nome);
  sprintf(a2, "el_%s_pfiso_myphoton03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_myphoton03, a2);

  sprintf(a1, "el_%s_pfiso_myneutral04", nome);
  sprintf(a2, "el_%s_pfiso_myneutral04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_myneutral04, a2);

  sprintf(a1, "el_%s_pfiso_mycharged04", nome);
  sprintf(a2, "el_%s_pfiso_mycharged04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_mycharged04, a2);
  
  sprintf(a1, "el_%s_pfiso_myphoton04", nome);
  sprintf(a2, "el_%s_pfiso_myphoton04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_myphoton04, a2);

  sprintf(a1, "el_%s_pfiso_neutral", nome);
  sprintf(a2, "el_%s_pfiso_neutral[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_neutral, a2);

  sprintf(a1, "el_%s_pfiso_charged", nome);
  sprintf(a2, "el_%s_pfiso_charged[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_charged, a2);
  
  sprintf(a1, "el_%s_pfiso_photon", nome);
  sprintf(a2, "el_%s_pfiso_photon[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_pfiso_photon, a2);

  sprintf(a1, "el_%s_hcaliso03", nome);
  sprintf(a2, "el_%s_hcaliso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hcaliso03, a2);

  sprintf(a1, "el_%s_ecaliso03", nome);
  sprintf(a2, "el_%s_ecaliso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ecaliso03, a2);
  
  sprintf(a1, "el_%s_tkiso03", nome);
  sprintf(a2, "el_%s_tkiso03[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_tkiso03, a2);

  sprintf(a1, "el_%s_hcaliso04", nome);
  sprintf(a2, "el_%s_hcaliso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hcaliso04, a2);

  sprintf(a1, "el_%s_hcalsolidiso04", nome);
  sprintf(a2, "el_%s_hcalsolidiso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_hcalsolidiso04, a2);

  //sprintf(a1, "el_%s_hcalbciso03", nome);
  //sprintf(a2, "el_%s_hcalbciso03[el_%s_n]/F", nome, nome);
  //tree->Branch(a1, &el_hcalbciso03, a2);

  //sprintf(a1, "el_%s_hcalbciso04", nome);
  //sprintf(a2, "el_%s_hcalbciso04[el_%s_n]/F", nome, nome);
  //tree->Branch(a1, &el_hcalbciso04, a2);

  sprintf(a1, "el_%s_ecaliso04", nome);
  sprintf(a2, "el_%s_ecaliso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ecaliso04, a2);
  
  sprintf(a1, "el_%s_tkiso04", nome);
  sprintf(a2, "el_%s_tkiso04[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_tkiso04, a2);
  
  sprintf(a1, "el_%s_tkdrv", nome);
  sprintf(a2, "el_%s_tkdrv[el_%s_n]/O", nome, nome);
  tree->Branch(a1, &el_tkdrv, a2);
 
  sprintf(a1, "el_%s_ecaldrv", nome);
  sprintf(a2, "el_%s_ecaldrv[el_%s_n]/O", nome, nome);
  tree->Branch(a1, &el_ecaldrv, a2);

  sprintf(a1, "el_%s_ip_ctf", nome);
  sprintf(a2, "el_%s_ip_ctf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ip_ctf, a2);

  sprintf(a1, "el_%s_ip_gsf", nome);
  sprintf(a2, "el_%s_ip_gsf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_ip_gsf, a2);

  sprintf(a1, "el_%s_dist", nome);
  sprintf(a2, "el_%s_dist[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dist, a2);

  sprintf(a1, "el_%s_dcot", nome);
  sprintf(a2, "el_%s_dcot[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_dcot, a2);

  sprintf(a1, "el_%s_hp_1pxb", nome);
  sprintf(a2, "el_%s_hp_1pxb[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_1pxb, a2);

  sprintf(a1, "el_%s_hp_1pxf", nome);
  sprintf(a2, "el_%s_hp_1pxf[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_1pxf, a2);

  sprintf(a1, "el_%s_conv", nome);
  sprintf(a2, "el_%s_conv[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_conv, a2);  

  sprintf(a1, "el_%s_regr_energy", nome);
  sprintf(a2, "el_%s_regr_energy[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_regr_energy, a2);

  sprintf(a1, "el_%s_regr_energyerr", nome);
  sprintf(a2, "el_%s_regr_energyerr[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_regr_energyerr, a2);

  sprintf(a1, "el_%s_eleopout", nome);
  sprintf(a2, "el_%s_eleopout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eleopout, a2);

  sprintf(a1, "el_%s_detaeleout", nome);
  sprintf(a2, "el_%s_detaeleout[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_detaeleout, a2);

  sprintf(a1, "el_%s_kfhits", nome);
  sprintf(a2, "el_%s_kfhits[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_kfhits, a2);

  sprintf(a1, "el_%s_kfchi2", nome);
  sprintf(a2, "el_%s_kfchi2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_kfchi2, a2);

  sprintf(a1, "el_%s_psenergy", nome);
  sprintf(a2, "el_%s_psenergy[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_psenergy, a2);  

  sprintf(a1, "el_%s_passmvapresel", nome);
  sprintf(a2, "el_%s_passmvapresel[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_passmvapresel, a2);

  sprintf(a1, "el_%s_passcutpresel", nome);
  sprintf(a2, "el_%s_passcutpresel[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_passcutpresel, a2);

  sprintf(a1, "el_%s_psenergypf", nome);
  sprintf(a2, "el_%s_psenergypf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_psenergypf, a2);  
  
  sprintf(a1, "el_%s_nbrempf", nome);
  sprintf(a2, "el_%s_nbrempf[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_nbrempf, a2);

  sprintf(a1, "el_%s_eseedpf", nome);
  sprintf(a2, "el_%s_eseedpf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_eseedpf, a2);

  sprintf(a1, "el_%s_epf", nome);
  sprintf(a2, "el_%s_epf[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_epf, a2);

  sprintf(a1, "el_%s_psly1", nome);
  sprintf(a2, "el_%s_psly1[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_psly1, a2);

  sprintf(a1, "el_%s_psly2", nome);
  sprintf(a2, "el_%s_psly2[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_psly2, a2);

  sprintf(a1, "el_%s_psnstriply1", nome);
  sprintf(a2, "el_%s_psnstriply1[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_psnstriply1, a2);

  sprintf(a1, "el_%s_psnstriply2", nome);
  sprintf(a2, "el_%s_psnstriply2[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_psnstriply2, a2);

  sprintf(a1, "el_%s_D0Vtx", nome);
  sprintf(a2, "el_%s_D0Vtx[el_%s_n][100]/F", nome, nome);
  tree->Branch(a1, &el_D0Vtx, a2);

  sprintf(a1, "el_%s_DZVtx", nome);
  sprintf(a2, "el_%s_DZVtx[el_%s_n][100]/F", nome, nome);
  tree->Branch(a1, &el_DZVtx, a2);

  sprintf(a1, "el_%s_must", nome);
  sprintf(a2, "el_%s_must[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_must, a2);

  sprintf(a1, "el_%s_mustnc", nome);
  sprintf(a2, "el_%s_mustnc[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_mustnc, a2);

  sprintf(a1, "el_%s_r9", nome);
  sprintf(a2, "el_%s_r9[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_r9, a2);

  sprintf(a1, "el_%s_gsfchi2", nome);
  sprintf(a2, "el_%s_gsfchi2[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_gsfchi2, a2);
  
  sprintf(a1, "el_%s_ip3d", nome);
  sprintf(a2, "el_%s_ip3d[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ip3d, a2);

  sprintf(a1, "el_%s_ip3d_err", nome);
  sprintf(a2, "el_%s_ip3d_err[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ip3d_err, a2);

  sprintf(a1, "el_%s_ip3d_sig", nome);
  sprintf(a2, "el_%s_ip3d_sig[el_%s_n]/I", nome, nome);
  tree->Branch(a1, &el_ip3d_sig, a2);

  //sprintf(a1, "el_%s_sc_time", nome);
  //sprintf(a2, "el_%s_sc_time[el_%s_n]/I", nome, nome);
  //tree->Branch(a1, &el_sc_time, a2);

  sprintf(a1, "el_%s_conv_vtxProb", nome);
  sprintf(a2, "el_%s_conv_vtxProb[el_%s_n]/F", nome, nome);
  tree->Branch(a1, &el_conv_vtxProb, a2);
}


bool GlobeElectrons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(electronColl, elH);

  edm::Handle<reco::GsfElectronCollection> calibEleH;
  edm::Handle<edm::ValueMap<double> > calibEnergyH;
  edm::Handle<edm::ValueMap<double> > calibEnergyErrH;
  edm::Handle<edm::ValueMap<double> > corrEnergyH;
  edm::Handle<edm::ValueMap<double> > corrEnergyErrH;
  edm::ValueMap<double>* calibEnergy=0;
  edm::ValueMap<double>* calibEnergyErr=0; 
  edm::ValueMap<double>* corrEnergy=0; 
  edm::ValueMap<double>* corrEnergyErr=0;

  if (applyEnergyCorrection) {
    iEvent.getByLabel("calibratedElectrons", "calibratedGsfElectrons", calibEleH);
    iEvent.getByLabel("calibratedElectrons" ,"eneRegForGsfEle", calibEnergyH);
    calibEnergy = const_cast<edm::ValueMap<double>* > (calibEnergyH.product());
    iEvent.getByLabel("calibratedElectrons" ,"eneErrorRegForGsfEle", calibEnergyErrH);
    calibEnergyErr = const_cast<edm::ValueMap<double>* > (calibEnergyErrH.product());
    iEvent.getByLabel("eleRegressionEnergy" ,"eneRegForGsfEle", corrEnergyH);
    corrEnergy = const_cast<edm::ValueMap<double>* > (corrEnergyH.product());
    iEvent.getByLabel("eleRegressionEnergy" ,"eneErrorRegForGsfEle", corrEnergyErrH);
    corrEnergyErr = const_cast<edm::ValueMap<double>* > (corrEnergyErrH.product());
  }

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(conversionColl, hConversions);

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamSpotColl, bsHandle);
  const reco::BeamSpot &thebs = *bsHandle.product();

  edm::Handle<reco::TrackCollection> tkH;
  iEvent.getByLabel(trackColl, tkH);

  edm::Handle<reco::GsfTrackCollection> tkH2;
  iEvent.getByLabel(trackColl2, tkH2);
  
  edm::Handle<reco::SuperClusterCollection> superClustersBarrelH; 
  iEvent.getByLabel(hybridSuperClusterColl,superClustersBarrelH);
  
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vertexColl, vtxH);


  //edm::Handle<double> rhoHandle;
  //iEvent.getByLabel(rhoCollection, rhoHandle);
  //double rho = *(rhoHandle.product());
  
  // transient track builder needed for ele ID MVA
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_);  
  const TransientTrackBuilder thebuilder = *(trackBuilder_.product());
  
  /*
  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel(pfColl, pfHandle);
  edm::Handle<reco::PFCandidateCollection> pfHandlePu;
  iEvent.getByLabel("pfPileUp", pfHandlePu);
  */

  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();

  edm::Handle< HBHERecHitCollection > hbhe ;
  iEvent.getByLabel(hcalHitColl, hbhe);

  edm::Handle<EERecHitCollection> pEERecHitH;
  edm::Handle<EBRecHitCollection> pEBRecHitH;
  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  const EcalRecHitCollection *barrelRecHits = pEBRecHitH.product();
  const EcalRecHitCollection *endcapRecHits = pEERecHitH.product();
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;

  const CaloSubdetectorGeometry *geometryES = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  CaloSubdetectorTopology *topology_p = 0;
  if (geometryES) topology_p = new EcalPreshowerTopology(geoHandle);
  
  ecalLazyTool = new EcalClusterLazyTools(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl);

  edm::Handle<EcalRecHitCollection> ESRecHits;
  iEvent.getByLabel(ecalHitESColl , ESRecHits);

  rechits_map_.clear();
  if (ESRecHits.isValid()) {
    EcalRecHitCollection::const_iterator it;
    for (it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
      // remove bad ES rechits
      if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
      //Make the map of DetID, EcalRecHit pairs
      rechits_map_.insert(std::make_pair(it->id(), *it));
    }
  }

  el_sc->Clear();
  el_p4->Clear(); 
  el_momvtx->Clear();
  el_momvtxconst->Clear();
  el_momcalo->Clear();
  el_posvtx->Clear();
  el_poscalo->Clear(); 
  el_catbased ->clear();  
  el_schits->clear();
  el_bchits->clear();
  el_n = 0;
  
  if (debug_level > 9)
    std::cout << "GlobeElectrons: Electron collection size: "<< elH->size() << std::endl;
  
  for(reco::GsfElectronCollection::const_iterator igsf = elH->begin(); igsf != elH->end(); igsf++) {

    if (el_n >= MAX_ELECTRONS) {
      std::cout << "GlobeElectrons: WARNING TOO MANY ELECTRONS: " << elH->size() << " (allowed " << MAX_ELECTRONS << ")" << std::endl;
      break;
    }
    
    reco::GsfElectron egsf = reco::GsfElectron(*igsf);
    
    if(gCUT->cut(egsf))
      continue;
    
    float phi = egsf.superCluster()->phi();
    float theta = (2*atan(exp(-egsf.superCluster()->eta())));
    float en = egsf.superCluster()->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);

    new ((*el_sc)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_sc->At(el_n))->SetXYZT(px, py, pz, en);
    
    new ((*el_p4)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_p4->At(el_n))->SetXYZT(egsf.px(), egsf.py(), egsf.pz(), egsf.energy());

    new ((*el_momvtxconst)[el_n]) TVector3();
    ((TVector3 *)el_momvtxconst->At(el_n))->SetXYZ(egsf.trackMomentumAtVtxWithConstraint().x(), 
                                              egsf.trackMomentumAtVtxWithConstraint().y(), egsf.trackMomentumAtVtxWithConstraint().z());
    
    new ((*el_momvtx)[el_n]) TVector3();
    ((TVector3 *)el_momvtx->At(el_n))->SetXYZ(egsf.trackMomentumAtVtx().x(), 
                                              egsf.trackMomentumAtVtx().y(), egsf.trackMomentumAtVtx().z());
    
    new ((*el_momcalo)[el_n]) TVector3();
    ((TVector3 *)el_momcalo->At(el_n))->SetXYZ(egsf.trackMomentumAtCalo().x(), 
                                               egsf.trackMomentumAtCalo().y(), egsf.trackMomentumAtCalo().z());

    new ((*el_momout)[el_n]) TVector3();
    ((TVector3 *)el_momout->At(el_n))->SetXYZ(egsf.trackMomentumOut().x(), 
                                               egsf.trackMomentumOut().y(), egsf.trackMomentumOut().z());
    
    new ((*el_posvtx)[el_n]) TVector3();
    ((TVector3 *)el_posvtx->At(el_n))->SetXYZ(egsf.trackPositionAtVtx().x(), 
                                              egsf.trackPositionAtVtx().y(), egsf.trackPositionAtVtx().z());

    new ((*el_poscalo)[el_n]) TVector3();
    ((TVector3 *)el_poscalo->At(el_n))->SetXYZ(egsf.trackPositionAtCalo().x(), 
                                               egsf.trackPositionAtCalo().y(), egsf.trackPositionAtCalo().z());
    
    el_pout[el_n] = egsf.trackMomentumOut().R();
    el_pin[el_n] = egsf.trackMomentumAtVtx().R();
    
    el_1pxb[el_n] = 0;
    el_1pxf[el_n] = 0;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel())
      el_1pxb[el_n] = 1;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelEndcap())
      el_1pxf[el_n] = 1;

    el_fbrem[el_n] = egsf.fbrem();

    // FIXME
    bool passconversionveto = !ConversionTools::hasMatchedConversion(egsf, hConversions, thebs.position());
    el_conv[el_n] = int(passconversionveto);

    el_eseedopout[el_n] = egsf.eSeedClusterOverPout();
    el_eseedopin[el_n] = egsf.eSeedClusterOverP();
    el_eopin[el_n] = egsf.eSuperClusterOverP();
    el_eleopout[el_n] = egsf.eEleClusterOverPout();
    el_detain[el_n] = egsf.deltaEtaSuperClusterTrackAtVtx();
    el_dphiin[el_n] = egsf.deltaPhiSuperClusterTrackAtVtx();
    el_detaout[el_n] = egsf.deltaEtaSeedClusterTrackAtCalo();
    el_dphiout[el_n] = egsf.deltaPhiSeedClusterTrackAtCalo();
    el_detaeleout[el_n] = egsf.deltaEtaEleClusterTrackAtCalo();
    

    el_nambtk[el_n] = egsf.ambiguousGsfTracksSize();
    el_class[el_n] = egsf.classification();

    el_nbrem[el_n] = egsf.numberOfBrems();
   
    el_e5x5[el_n] = egsf.e5x5();
    el_e2x5[el_n] = egsf.e2x5Max();
    el_e1x5[el_n] = egsf.e1x5();

    el_sieie[el_n] = egsf.sigmaIetaIeta();

    //el_1oe_1op[el_n] = 1./egsf.ecalEnergy() - 1./egsf.trackMomentumAtVtx().R();
    el_1oe_1op[el_n] = 1./egsf.ecalEnergy() - 1./egsf.p();
    el_r9[el_n] = ecalLazyTool.e3x3(*(egsf.superCluster()->seed())) / egsf.superCluster()->rawEnergy();

    //el_sc_time[el_n] = ecalLazyTool.SuperClusterTime(*(egsf.superCluster()), iEvent);

    el_must[el_n] = -9999.;
    el_mustnc[el_n] = -1;
    reco::Mustache m;
    m.MustacheID(*(egsf.superCluster()), el_mustnc[el_n], el_must[el_n]);

    // Regression Correction
    if (!ecorr_.IsInitialized()) {
      if (energyCorrectionsFromDB and regressionVersion == "V3") {
	std::cout << "DB version available only for V2" << std::endl;
	energyCorrectionsFromDB = false;
      }
      
      if (!energyCorrectionsFromDB) {
	char filename[500];
	char* descr = getenv("CMSSW_BASE");
	sprintf(filename, "%s/src/HiggsAnalysis/HiggsTo2photons/data/%s", descr, energyRegFilename.c_str());
	if (fexist(filename)) {
	  sprintf(filename, "http://cern.ch/sani/%s", energyRegFilename.c_str());
	  ecorr_.Initialize(iSetup, filename);
	} else {
	  ecorr_.Initialize(iSetup, filename);
	} 
      } else {
	ecorr_.Initialize(iSetup, "wgbrph", true); //FIXME no wgbrele in DB
      }
    }

    //if (regressionVersion == "V3") {
    //std::pair<double,double> cor = ecorr_.CorrectedEnergyWithErrorV3(egsf, *vtxH, rho, ecalLazyTool, iSetup);
    //el_regr_energy[el_n]    = cor.first;
    //el_regr_energyerr[el_n] = cor.second;
    //} else {
    std::pair<double,double> cor = ecorr_.CorrectedEnergyWithError(egsf, ecalLazyTool);
    el_regr_energy[el_n]    = cor.first;
    el_regr_energyerr[el_n] = cor.second;
    //}
    
    reco::GsfElectronRef myElectronRef(elH, std::distance(elH->begin(), igsf));
    if (applyEnergyCorrection) {
      reco::GsfElectronRef calibEleRef(calibEleH, std::distance(elH->begin(), igsf));
      el_corr_energy[el_n]     = (*corrEnergy)[myElectronRef];
      el_corr_energyerr[el_n]  = (*corrEnergyErr)[myElectronRef];
      el_calib_energy[el_n]    = (*calibEnergy)[calibEleRef];
      el_calib_energyerr[el_n] = (*calibEnergyErr)[calibEleRef];
      
      new ((*el_p4_corr)[el_n]) TLorentzVector();
      ((TLorentzVector *)el_p4_corr->At(el_n))->SetXYZT(calibEleRef->px(), calibEleRef->py(), calibEleRef->pz(), calibEleRef->energy());
    }

    // ES variables
    el_eseffsixix[el_n] = 0.;
    el_eseffsiyiy[el_n] = 0.;
    if (ESRecHits.isValid() && (fabs(egsf.superCluster()->eta()) > 1.6 && fabs(egsf.superCluster()->eta()) < 3)) {
      std::vector<float> elESHits0 = gES->getESHits(egsf.superCluster()->x(), egsf.superCluster()->y(), egsf.superCluster()->z(), rechits_map_, geometry, topology_p, 0);
      std::vector<float> elESShape = gES->getESShape(elESHits0);
      el_eseffsixix[el_n] = elESShape[0];
      el_eseffsiyiy[el_n] = elESShape[1];
    }

    el_hoe[el_n] = egsf.hcalOverEcal();
    el_hoed1[el_n] = egsf.hcalDepth1OverEcal();
    el_hoed2[el_n] = egsf.hcalDepth2OverEcal();

    //el_hoebc[el_n] = egsf.hcalOverEcalBc();
    //el_hoebcd1[el_n] = egsf.hcalDepth1OverEcalBc();
    //el_hoebcd2[el_n] = egsf.hcalDepth2OverEcalBc();

    el_d0[el_n] = egsf.gsfTrack()->d0();
    el_z0[el_n] = egsf.gsfTrack()->dz(); 
    el_chi2[el_n] = egsf.gsfTrack()->chi2();
    el_dof[el_n] = egsf.gsfTrack()->ndof();
    el_validhits[el_n] = egsf.gsfTrack()->numberOfValidHits();
    el_losthits[el_n] = egsf.gsfTrack()->numberOfLostHits();

    math::XYZPoint vtxPoint(0.0,0.0,0.0);
    if (vtxH->size() != 0) {
      reco::VertexRef vtx(vtxH, 0);
      vtxPoint = math::XYZPoint(vtx->x(),vtx->y(),vtx->z());
    }

    el_ip_gsf[el_n] = egsf.gsfTrack()->dxy(vtxPoint);

    if (!doAodSim && (trackColl2.encode() != "electronGsfTracks")) {
      for(unsigned int j=0; j<tkH2->size(); j++) { 
        reco::GsfTrackRef tk2(tkH2, j);
        std::pair<unsigned int, float> result = sharedHits(*tk2, *(egsf.gsfTrack()));
        
        if (result.second > .999) {
          el_hp_expin[el_n] = tk2->trackerExpectedHitsInner().numberOfHits();
          el_hp_expout[el_n] = tk2->trackerExpectedHitsOuter().numberOfHits();
          break;
        }
      }
    } else {
      el_hp_expin[el_n] = egsf.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      el_hp_expout[el_n] = egsf.gsfTrack()->trackerExpectedHitsOuter().numberOfHits();
    }
 
    el_charge[el_n] = egsf.charge();
    
    el_ch_scpix[el_n] = egsf.scPixCharge();
    el_ch_gsf[el_n] = egsf.gsfTrack()->charge();

    if (egsf.isEB())
      el_crack[el_n] = 0;

    if (egsf.isEE())
      el_crack[el_n] = 1;

    if (egsf.isEBEEGap())
      el_crack[el_n] = 3;

    if (egsf.isEBGap())
      el_crack[el_n] = 4;

    if (egsf.isEEGap())
      el_crack[el_n] = 5;

    //int cmsTkind = -1;

    el_tkind[el_n] = -1;
    el_scind[el_n] = -1;
  
    if(egsf.closestCtfTrackRef().isNonnull()) {
      el_kfhits[el_n] = egsf.closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement();
      el_kfchi2[el_n] = egsf.closestCtfTrackRef()->normalizedChi2();
      for(unsigned int j=0; j<tkH->size(); ++j) {
        reco::TrackRef tk(tkH, j);
        if(gCUT->cut(*tk))
          continue; 
        if (tk == egsf.closestCtfTrackRef()) {
          el_tkind[el_n] = j;
          el_ip_ctf[el_n] = egsf.closestCtfTrackRef()->dxy(vtxPoint);
        }
      }
    } else {
      el_kfhits[el_n] = -1;
      el_kfchi2[el_n] = -999;
      el_tkind[el_n] = -1;
      el_ip_ctf[el_n] = 0;
    }
    
    int index = 0;
    // loop over the two SC collections
    for(int z = 0; z<2; ++z) {
      if (z == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersBarrelH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersBarrelH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))
	    continue;
          // passed cuts

          if (&(*egsf.superCluster()) == &(*cluster)) {
            el_scind[el_n] = index; 
            el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*barrelRecHits), &(*topology))[0]);
            std::vector<float> vCov = ecalLazyTool->localCovariances(*(cluster->seed()));
            el_sipip[el_n] = sqrt(vCov[2]);
            break;
          }
          index++;
        }
      }
      
      if (z == 1) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersEndcapH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))continue;
          // passed cuts 

          if (&(*(egsf.superCluster())) == &(*cluster)) {
            el_scind[el_n] = index;
            el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*endcapRecHits), &(*topology))[0]);
            std::vector<float> vCov = ecalLazyTool->localCovariances(*(cluster->seed()));
            el_sipip[el_n] = sqrt(vCov[2]);
            break;
          }
          index++;
        }
      }
    }

    if (el_scind[el_n] == -1) {
      if(fabs(egsf.superCluster()->eta()) < 1.479) {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*barrelRecHits), &(*topology))[0]);
      } else {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*endcapRecHits), &(*topology))[0]);
      }
    }

    std::vector<UInt_t> schits, bchits;    
    for (reco::CaloCluster_iterator clusterIt=egsf.superCluster()->clustersBegin(); clusterIt!=egsf.superCluster()->clustersEnd(); clusterIt++) {
      for (unsigned int j=0; j<(*clusterIt)->hitsAndFractions().size(); j++)
	schits.push_back(((*clusterIt)->hitsAndFractions())[j].first);
      el_schits->push_back(schits);
    }
    
    for (unsigned int j=0; j<egsf.superCluster()->seed()->hitsAndFractions().size(); j++)
      bchits.push_back((egsf.superCluster()->seed()->hitsAndFractions())[j].first);
    el_bchits->push_back(bchits);

    //el_mva_noiso[el_n] = egsf.mva();
    //el_mva[el_n] = mvaEstimator->mva(egsf, vtxH->size());
    
    el_mva_nontrig[el_n] = myMVANonTrig->mvaValue(egsf, 
						    vtxH->front(), 
						    thebuilder,
						    ecalLazyTool,
						    false);
       
    el_mva_trig[el_n] = myMVATrig->mvaValue(egsf,
					      vtxH->front(),
					      thebuilder,
					      ecalLazyTool,
					      false);
    
    //std::cout<<"el_n el_mva_trig el_mva_nontrig "<<el_n<<" "<<el_mva_trig[el_n]<<" "<<el_mva_nontrig[el_n]<<std::endl;

    /* NEW variables */
    
    if (!egsf.pflowSuperCluster().isNull()) {
      el_nbrempf[el_n] = egsf.pflowSuperCluster()->clustersSize();
      el_eseedpf[el_n] = egsf.pflowSuperCluster()->seed()->energy();
      el_epf[el_n] = egsf.pflowSuperCluster()->energy();
      el_psenergypf[el_n] = egsf.pflowSuperCluster()->preshowerEnergy();
    } else {
      el_nbrempf[el_n] = -1;
      el_eseedpf[el_n] = -1;
      el_epf[el_n] = -1;
      el_psenergypf[el_n] = -1;
    }
    
    el_passcutpresel[el_n] = (Int_t)egsf.passingCutBasedPreselection();
    el_passmvapresel[el_n] = (Int_t)egsf.passingMvaPreselection();
    el_psenergy[el_n] = egsf.superCluster()->preshowerEnergy();  

    el_psly1[el_n] = 0;
    el_psly2[el_n] = 0;
    el_psnstriply1[el_n] = 0;
    el_psnstriply2[el_n] = 0;
    for (reco::CaloCluster_iterator it=egsf.superCluster()->preshowerClustersBegin();
	 it != egsf.superCluster()->preshowerClustersEnd(); it++) {

      if (ESDetId(((*it)->hitsAndFractions())[0].first).plane() == 1) {
	el_psly1[el_n] += (*it)->energy();
	el_psnstriply1[el_n] += (*it)->hitsAndFractions().size();
      } else {
	el_psly2[el_n] += (*it)->energy();
	el_psnstriply2[el_n] += (*it)->hitsAndFractions().size();
      }
    }

    el_ecaldrv[el_n] = egsf.ecalDrivenSeed();
    el_tkdrv[el_n] = egsf.trackerDrivenSeed();
    
    el_pfiso_charged[el_n] = egsf.pfIsolationVariables().chargedHadronIso;
    el_pfiso_photon[el_n] = egsf.pfIsolationVariables().photonIso;
    el_pfiso_neutral[el_n] = egsf.pfIsolationVariables().neutralHadronIso;

    el_tkiso04[el_n] = egsf.dr04TkSumPt();
    el_ecaliso04[el_n] = egsf.dr04EcalRecHitSumEt();
    el_hcaliso04[el_n] = egsf.dr04HcalTowerSumEt();    

    el_tkiso03[el_n] = egsf.dr03TkSumPt();
    el_ecaliso03[el_n] = egsf.dr03EcalRecHitSumEt();
    el_hcaliso03[el_n] = egsf.dr03HcalTowerSumEt();

    EgammaTowerIsolation hcaliso03(0.3, 0., 0, -1, ctH.product());
    EgammaTowerIsolation hcaliso04(0.4, 0., 0, -1, ctH.product());

    el_hcalsolidiso03[el_n] = hcaliso03.getTowerEtSum(&(*(egsf.superCluster())) );
    el_hcalsolidiso04[el_n] = hcaliso04.getTowerEtSum(&(*(egsf.superCluster())) );

    //el_hcalbciso03[el_n] = egsf.dr03HcalTowerSumEtBc();
    //el_hcalbciso04[el_n] = egsf.dr04HcalTowerSumEtBc();

    // Fill out electron identification
    std::vector<edm::Handle<edm::ValueMap<float> > > eIDVM(9); 
    std::vector<int> results;
    
    for(unsigned int j=0; j<eIDLabels.size(); j++) {
      if (iEvent.getByLabel(eIDLabels[j], eIDVM[j])) {  
        const edm::ValueMap<float>& eIDmapTemp = *eIDVM[j];
        reco::GsfElectronRef electronRef(elH, std::distance(elH->begin(), igsf));
        results.push_back((Int_t)eIDmapTemp[electronRef]);
      } else {
        results.push_back(-1);
      }
    }
    
    el_catbased->push_back(results);

    el_dist[el_n] = (egsf.convDist() == -9999.? 9999:egsf.convDist());
    el_dcot[el_n] = (egsf.convDcot() == -9999.? 9999:egsf.convDcot());
  
    // loop through vertices for d0 and dZ w.r.t. each vertex
    // need number of vertices and vertices' positions
    std::cout<<"now I want d0 dZ:  eln "<<el_n<<std::endl;
    int maxV = std::max(100, (int)vtxH->size());
    std::cout<<"maxV (int)vtxH->size() "<<maxV<<"  "<<(int)vtxH->size()<<std::endl;
    for(int iv=0; iv<maxV; iv++){
      std::cout<<"gettting vertexref:  iv "<<iv<<std::endl;
      reco::VertexRef v(vtxH, iv);
      math::XYZPoint vtxPoint = math::XYZPoint(v->x(), v->y(), v->z());
      std::cout<<"got vertexref and vtxPoint:  iv "<<iv<<std::endl;
      
      el_D0Vtx[el_n][iv] = egsf.gsfTrack()->dxy(vtxPoint);
      std::cout<<"el_D0Vtx["<<el_n<<"]["<<iv<<"]  "<<el_D0Vtx[el_n][iv]<<std::endl;
      el_DZVtx[el_n][iv] = egsf.gsfTrack()->dz(vtxPoint);
      std::cout<<"el_DZVtx["<<el_n<<"]["<<iv<<"]  "<<el_DZVtx[el_n][iv]<<std::endl;
    }

    el_n++;
  }
  
  return true;
}

std::pair<unsigned int, float> GlobeElectrons::sharedHits(const reco::Track& trackA, const reco::Track& trackB) {
  
  unsigned int shared = 0;
  for(trackingRecHit_iterator tkHitA = trackA.recHitsBegin(); tkHitA !=trackA.recHitsEnd(); ++tkHitA){
    for(trackingRecHit_iterator tkHitB = trackB.recHitsBegin();
        tkHitB !=trackB.recHitsEnd(); ++tkHitB){
      if( (**tkHitA).isValid() && (**tkHitB).isValid() &&(**tkHitA).sharesInput( &(**tkHitB),TrackingRecHit::all)) {
        shared++;
        break;
      }
    }
  }

  float fraction = (float) shared/std::min(trackA.found(),trackB.found());
  return std::make_pair(shared,fraction);
}








