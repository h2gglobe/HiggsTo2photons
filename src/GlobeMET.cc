#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMET.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

GlobeMET::GlobeMET(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  caloMETColl =  iConfig.getParameter<edm::InputTag>("CaloMETColl");
  if(strcmp(nome, "tcmet") == 0) 
    tcMETColl =  iConfig.getParameter<edm::InputTag>("TcMETColl");
  
  if(strcmp(nome, "pfmet") == 0) {
    pfMETColl =  iConfig.getParameter<edm::InputTag>("PFMETColl");
    pfType1METColl = iConfig.getParameter<edm::InputTag>("PFMETTYPE1Coll");
  }

  muonGlobalColl =  iConfig.getParameter<edm::InputTag>("MuonColl");

  jetColl =  iConfig.getParameter<edm::InputTag>("JetColl_algo1");

  debug_level = iConfig.getParameter<int>("Debug_Level");
}

void GlobeMET::defineBranch(TTree* tree) {
   
  if(strcmp(nome, "tcmet") != 0 && strcmp(nome, "pfmet") != 0){
    tree->Branch("met_met", &met_met, "met_met/F");
    tree->Branch("met_phi", &met_phi, "met_phi/F");
    tree->Branch("met_met_nocalo", &met_met_nocalo, "met_met_nocalo/F");
    tree->Branch("met_phi_nocalo", &met_phi_nocalo, "met_phi_nocalo/F");
    tree->Branch("met_met_crossed", &met_met_crossed, "met_met_crossed/F");
    tree->Branch("met_phi_crossed", &met_phi_crossed, "met_phi_crossed/F");
    tree->Branch("met_met_s9", &met_met_s9, "met_met_s9/F");
    tree->Branch("met_phi_s9", &met_phi_s9, "met_phi_s9/F");
    tree->Branch("met_met_mip", &met_met_mip, "met_met_mip/F");
    tree->Branch("met_phi_mip", &met_phi_mip, "met_phi_mip/F");
    tree->Branch("met_met_jet", &met_met_jet, "met_met_jet/F");
    tree->Branch("met_phi_jet", &met_phi_jet, "met_phi_jet/F");
  }

  if(strcmp(nome, "tcmet") == 0) {
    tree->Branch("met_tcmet", &met_tcmet, "met_tcmet/F");
    tree->Branch("met_phi_tcmet", &met_phi_tcmet, "met_phi_tcmet/F");
  }

  if(strcmp(nome, "pfmet") == 0) {
    tree->Branch("met_pfmet", &met_pfmet, "met_pfmet/F");
    tree->Branch("met_phi_pfmet", &met_phi_pfmet, "met_phi_pfmet/F");
    tree->Branch("met_sumet_pfmet", &met_sumet_pfmet, "met_sumet_pfmet/F");
    tree->Branch("met_mEtSig_pfmet", &met_mEtSig_pfmet, "met_mEtSig_pfmet/F");
    tree->Branch("met_significance_pfmet", &met_significance_pfmet, "met_significance_pfmet/F");
    
    tree->Branch("met_pfmetType1", &met_pfmetType1, "met_pfmetType1/F");
    tree->Branch("met_phi_pfmetType1", &met_phi_pfmetType1, "met_phi_pfmetType1/F");
    tree->Branch("met_sumet_pfmetType1", &met_sumet_pfmetType1, "met_sumet_pfmetType1/F");
    tree->Branch("met_mEtSig_pfmetType1", &met_mEtSig_pfmetType1, "met_mEtSig_pfmetType1/F");
    tree->Branch("met_significance_pfmetType1", &met_significance_pfmetType1, "met_significance_pfmetType1/F");
  }
}

bool GlobeMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if(strcmp(nome, "tcmet") != 0 && strcmp(nome, "pfmet") != 0){
    if(debug_level > 99 ) 
      std::cout << "GlobeMET: Start " << std::endl;
    
    edm::Handle<reco::CaloMETCollection> metH;
    iEvent.getByLabel(caloMETColl, metH);
   
    met_met = metH->begin()->et();
    met_phi = metH->begin()->phi();
    
    met_met_nocalo = met_met;    
    met_phi_nocalo = met_phi;
    correctMETmuons(iEvent, met_met_nocalo, met_phi_nocalo, NoCaloCorrection);
    
    met_met_crossed = met_met;
    met_phi_crossed = met_phi;
    correctMETmuons(iEvent, met_met_crossed, met_phi_crossed, CrossedEnergyCorrection);
    
    met_met_s9 = met_met;
    met_phi_s9 = met_phi;
    correctMETmuons(iEvent, met_met_s9, met_phi_s9, S9EnergyCorrection);
    
    met_met_mip = met_met;
    met_phi_mip = met_phi;
    correctMETmuons(iEvent, met_met_mip, met_phi_mip, ExpectedMipEnergyCorrection);
    
    met_met_jet = met_met;
    met_phi_jet = met_phi;
    correctedJetMET(iEvent, met_met_jet, met_phi_jet, 5.);
    
    if (debug_level > 9)
      std::cout << "GlobeMet: MET collection size: "<< metH->size() << std::endl;
    
    if(debug_level > 99 ) 
      std::cout << "GlobeMET: End " << std::endl;
    
    return true;
  }

  if(strcmp(nome, "tcmet") == 0) {
    
    edm::Handle<reco::METCollection> tcmet_h;
    iEvent.getByLabel(tcMETColl, tcmet_h);
    
    met_tcmet = tcmet_h->begin()->et();
    met_phi_tcmet = tcmet_h->begin()->phi();

    return true;
  }

  if(strcmp(nome, "pfmet") == 0) {
    
    edm::Handle<reco::PFMETCollection> pfmet_h;
    iEvent.getByLabel(pfMETColl, pfmet_h);
    
    met_pfmet              = pfmet_h->begin()->et();
    met_phi_pfmet          = pfmet_h->begin()->phi();
    met_sumet_pfmet        = pfmet_h->begin()->sumEt();
    met_mEtSig_pfmet       = pfmet_h->begin()->mEtSig();
    met_significance_pfmet = pfmet_h->begin()->significance();
    
    iEvent.getByLabel(pfType1METColl, pfmet_h);
  
    met_pfmetType1              = pfmet_h->begin()->et();
    met_phi_pfmetType1          = pfmet_h->begin()->phi();
    met_sumet_pfmetType1        = pfmet_h->begin()->sumEt();
    met_mEtSig_pfmetType1       = pfmet_h->begin()->mEtSig();	      
    met_significance_pfmetType1 = pfmet_h->begin()->significance();

    //std::cout<<"met_pfmet  met_phi_pfmet  met_sumet_pfmet  met_mEtSig_pfmet  met_significance_pfmet  "<<
    //        met_pfmet<<"  "<<
    //        met_phi_pfmet<<"  "<<
    //        met_sumet_pfmet<<"  "<<
    //        met_mEtSig_pfmet<<"  "<<
    //        met_significance_pfmet<<std::endl;
    //
    //std::cout<<"met_pfmetType1  met_phi_pfmetType1  met_sumet_pfmetType1  met_mEtSig_pfmetType1  met_significance_pfmetType1  "<<
    //        met_pfmetType1<<"  "<<
    //        met_phi_pfmetType1<<"  "<<
    //        met_sumet_pfmetType1<<"  "<<
    //        met_mEtSig_pfmetType1<<"  "<<
    //        met_significance_pfmetType1<<std::endl;
    return true;
  }

  return true;

}
// correct MET energies for Muons
// (input parameters are corrected by the algorithm)
void GlobeMET::correctMETmuons(const edm::Event& iEvent, float& met, float& metPhi, CorrectionType type) {
  
  // first, account for muon momentum
  // take muon collection
  edm::Handle<reco::MuonCollection> metMuH;
  iEvent.getByLabel(muonGlobalColl, metMuH);
  
  double metx =  met*std::cos(metPhi);
  double mety =  met*std::sin(metPhi);
  for (reco::MuonCollection::const_iterator cand = metMuH->begin(); cand != metMuH->end(); ++cand ) {
    double pt0 = cand->pt(); 
    double phi0 = cand->phi(); 
    metx -= pt0*std::cos(phi0);
    mety -= pt0*std::sin(phi0);
  }
  met = std::sqrt(metx*metx+mety*mety);
  metPhi = std::atan2(mety, metx);
  
  if (type == NoCaloCorrection) 
    return;
  
  double muEx = 0.0;
  double muEy = 0.0;
  
  for (reco::MuonCollection::const_iterator mu = metMuH->begin(); mu != metMuH->end(); ++mu) {

    double theta = mu->theta();
    double phi = mu->phi();
    switch (type) {
    case CrossedEnergyCorrection:
      muEx += ( mu->calEnergy().em + mu->calEnergy().had + mu->calEnergy().ho )*sin(theta)*cos( phi );
      muEy += ( mu->calEnergy().em + mu->calEnergy().had + mu->calEnergy().ho )*sin(theta)*sin( phi );
      break;
    case S9EnergyCorrection:
      muEx += ( mu->calEnergy().emS9 + mu->calEnergy().hadS9 + mu->calEnergy().hoS9 )*sin(theta)*cos( phi );
      muEy += ( mu->calEnergy().emS9 + mu->calEnergy().hadS9 + mu->calEnergy().hoS9 )*sin(theta)*sin( phi );
      break;
    case ExpectedMipEnergyCorrection:
      // numbers are essential a wild guess
      if ( fabs(mu->eta()) < 1.5) { 
        // barrel
        muEx += (0.3 + 3.0 + 1.0)*sin(theta)*cos(phi);
        muEy += (0.3 + 3.0 + 1.0)*sin(theta)*sin(phi);
      } else {
        // endcap
        muEx += (0.35 + 3.5)*sin(theta)*cos(phi);
        muEy += (0.35 + 3.5)*sin(theta)*sin(phi);
      }
      break;
    default:
      std::cout << "Uknown MET correction type. Abort" << std::endl;
      assert(0);
    }
  }

  metx = met*cos(metPhi) + muEx;
  mety = met*sin(metPhi) + muEy;
  met   = std::sqrt(metx*metx + mety*mety);
  metPhi = atan2(mety, metx);
}

// list of jets must be supplied
void GlobeMET::correctedJetMET(const edm::Event& iEvent, float& met, float& metPhi, const float min_pt) {
  
  // first, account for muon momentum
  // take muon collection
  edm::Handle<reco::CaloJetCollection> jetH;
  iEvent.getByLabel(jetColl, jetH);

  //iterate over candidates, cast them to calojets and then correct for the energy
  double METX_uncorr = met*cos(metPhi);
  double METY_uncorr = met*sin(metPhi);  
  double Ex = 0.0;
  double Ey = 0.0;
  
  reco::CaloJetCollection::const_iterator cand_iter;
  for(reco::CaloJetCollection::const_iterator jet = jetH->begin(); jet != jetH->end(); ++jet) {
    int corr_factor = 1; // CHECK jet correction
    //jet correction doesn't do so well for recoJet pt < 30.0 
    if(jet->pt() > min_pt ) {  
      if(corr_factor > 0 ) { //args.jetcorrection is -999 in case of an error
        // std::cout << "Jet Pt: " << jet->pt() << "\t Jet phi: " << jet->phi() << "\tcorr: " << corr_factor <<std::endl;
        Ex = Ex + (corr_factor-1)*(jet->et())*cos(jet->phi());
        Ey = Ey + (corr_factor-1)*(jet->et())*sin(jet->phi());
      }
    }
  }

  double metx = METX_uncorr - Ex;
  double mety = METY_uncorr - Ey;
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
}
