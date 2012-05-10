#include <string>
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include <iostream>

GlobeJets::GlobeJets(const edm::ParameterSet& iConfig, const char* n = "algo1"): nome(n) {
  
  char a[100];
  sprintf (a,"JetColl_%s", nome);
  jetColl =  iConfig.getParameter<edm::InputTag>(a);
  calotowerColl =  iConfig.getParameter<edm::InputTag>("CaloTowerColl");
  trackColl =  iConfig.getParameter<edm::InputTag>("TrackColl");
  pfak5corrdata =  iConfig.getUntrackedParameter<std::string>("JetCorrectionData_"+std::string(n), "");
  pfak5corrmc   =  iConfig.getUntrackedParameter<std::string>("JetCorrectionMC_"+std::string(n), "");
  vertexColl = iConfig.getParameter<edm::InputTag>("VertexColl_std");
  jetMVAAlgos = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet> >("puJetIDAlgos_"+std::string(n), std::vector<edm::ParameterSet>(0));
  pfLooseId  = iConfig.getParameter<edm::ParameterSet>("pfLooseId");
  std::string strnome = nome;
    
  mvas_.resize(jetMVAAlgos.size());
  wp_levels_.resize(jetMVAAlgos.size());
  algos_.resize(jetMVAAlgos.size());

  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
      algos_[imva] = new PileupJetIdAlgo((jetMVAAlgos.at(imva)));
  }


  if (strnome.find("PF",0) == std::string::npos){
    sprintf (a,"JetTrackAssociationColl_%s", nome);
    jetTkAssColl =  iConfig.getParameter<edm::InputTag>(a);
  }

  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
  jet_tkind =  new std::vector<std::vector<unsigned short> >;
  jet_calotwind =  new std::vector<std::vector<unsigned short> >;
}

void GlobeJets::defineBranch(TTree* tree) {

  jet_p4 = new TClonesArray("TLorentzVector", MAX_JETS);
  
  char a1[50], a2[50];
  
  sprintf(a1, "jet_%s_n", nome);
  sprintf(a2, "jet_%s_n/I", nome);
  tree->Branch(a1, &jet_n, a2);
  
  sprintf(a1, "jet_%s_p4", nome);
  tree->Branch(a1, "TClonesArray", &jet_p4, 32000, 0);
  
  sprintf(a1, "jet_%s_emfrac", nome);
  sprintf(a2, "jet_%s_emfrac[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_emfrac, a2);

  sprintf(a1, "jet_%s_hadfrac", nome);
  sprintf(a2, "jet_%s_hadfrac[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_hadfrac, a2);

  sprintf(a1, "jet_%s_ntk", nome);
  sprintf(a2, "jet_%s_ntk[jet_%s_n]/I", nome, nome);
  tree->Branch(a1, &jet_ntk, a2);

  sprintf(a1, "jet_%s_erescale", nome);
  sprintf(a2, "jet_%s_erescale[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_erescale, a2);

  sprintf(a1, "jet_%s_dRMean", nome);
  sprintf(a2, "jet_%s_dRMean[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_dRMean, a2);

  sprintf(a1, "jet_%s_frac01", nome);
  sprintf(a2, "jet_%s_frac01[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac01, a2);

  sprintf(a1, "jet_%s_frac02", nome);
  sprintf(a2, "jet_%s_frac02[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac02, a2);

  sprintf(a1, "jet_%s_frac03", nome);
  sprintf(a2, "jet_%s_frac03[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac03, a2);

  sprintf(a1, "jet_%s_frac04", nome);
  sprintf(a2, "jet_%s_frac04[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac04, a2);

  sprintf(a1, "jet_%s_frac05", nome);
  sprintf(a2, "jet_%s_frac05[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac05, a2);

  sprintf(a1, "jet_%s_frac06", nome);
  sprintf(a2, "jet_%s_frac06[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac06, a2);

  sprintf(a1, "jet_%s_frac07", nome);
  sprintf(a2, "jet_%s_frac07[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_frac07, a2);

  sprintf(a1, "jet_%s_nNeutrals", nome);
  sprintf(a2, "jet_%s_nNeutrals[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_nNeutrals, a2);

  sprintf(a1, "jet_%s_beta", nome);
  sprintf(a2, "jet_%s_beta[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_beta, a2);

  sprintf(a1, "jet_%s_betaStar", nome);
  sprintf(a2, "jet_%s_betaStar[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_betaStar, a2);

  sprintf(a1, "jet_%s_dZ", nome);
  sprintf(a2, "jet_%s_dZ[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_dZ, a2);

  sprintf(a1, "jet_%s_nCharged", nome);
  sprintf(a2, "jet_%s_nCharged[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_nCharged, a2);

  sprintf(a1, "jet_%s_dR2Mean", nome);
  sprintf(a2, "jet_%s_dR2Mean[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_dR2Mean, a2);

  sprintf(a1, "jet_%s_betaStarClassic", nome);
  sprintf(a2, "jet_%s_betaStarClassic[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_betaStarClassic, a2);

  sprintf(a1, "jet_%s_pfloose", nome);
  sprintf(a2, "jet_%s_pfloose[jet_%s_n]/O", nome, nome);
  tree->Branch(a1, &jet_pfloose, a2);
 
  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
    
      std::string mvalabel = (jetMVAAlgos.at(imva)).getParameter<std::string>("label");

      mvas_[imva] = new float[MAX_JETS];
      wp_levels_[imva] = new int[MAX_JETS];

      sprintf(a1, "jet_%s_%s_mva", nome, mvalabel.c_str());
      sprintf(a2, "jet_%s_%s_mva[jet_%s_n]/F", nome, mvalabel.c_str(), nome);
      tree->Branch(a1, &mvas_[imva][0], a2);

      sprintf(a1, "jet_%s_%s_wp_level", nome, mvalabel.c_str());
      sprintf(a2, "jet_%s_%s_wp_level[jet_%s_n]/I", nome, mvalabel.c_str(), nome);
      tree->Branch(a1, &wp_levels_[imva][0], a2);
  }


  sprintf(a1, "jet_%s_pull_dy", nome);
  sprintf(a2, "jet_%s_pull_dy[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_pull_dy, a2);
  
  sprintf(a1, "jet_%s_pull_dphi", nome);
  sprintf(a2, "jet_%s_pull_dphi[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_pull_dphi, a2);
  
  sprintf(a1, "jet_%s_tkind", nome); 
  tree->Branch(a1, "std::vector<std::vector<unsigned short> >", &jet_tkind);

  sprintf(a1, "jet_%s_calotwind", nome);
  tree->Branch(a1, "std::vector<std::vector<unsigned short> >", &jet_calotwind);
}

bool GlobeJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::string strnome = nome;

  if (strnome.find("PF",0) == std::string::npos) {
    edm::Handle<reco::CaloJetCollection> jetH;
    iEvent.getByLabel(jetColl, jetH);  
    
    // take collections
        
    edm::Handle<CaloTowerCollection> ctH; 
    iEvent.getByLabel(calotowerColl, ctH);
    
    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(trackColl, tkH);
    
    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vertexColl, vtxH);
    
    jet_p4->Clear();
    jet_tkind->clear();
    jet_calotwind->clear();

    jet_n = 0;
    
    if (debug_level > 9)
      std::cout << "GlobeJets: Jet collection size: "<< jetH->size() << std::endl;
      
    // check if collection is present
    for(unsigned int i=0; i<jetH->size(); i++) {
      if (jet_n >= MAX_JETS) {
        std::cout << "GlobeJets: WARNING TOO MANY JETS: " << jetH->size() << " (allowed " << MAX_JETS << ")" << std::endl;
        break;
      }
      
      reco::CaloJetRef j(jetH, i);
	
      // apply the cuts
      if (gCUT->cut(*j)) continue;
      // passed cuts
      
      new ((*jet_p4)[jet_n]) TLorentzVector();
      ((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(j->px(), j->py(), j->pz(), j->energy()); 
      jet_emfrac[jet_n] = j->emEnergyFraction();
      jet_hadfrac[jet_n] = j->energyFractionHadronic();
      
      jet_erescale[jet_n] = 1;

      // Tracks and CaloTowers
      std::vector<CaloTowerPtr> towers = j->getCaloConstituents();
      std::vector<CaloTowerPtr>::const_iterator it;
      
      std::vector<unsigned short> temp;
      for(it = towers.begin(); it != towers.end(); ++it) {
        for(unsigned int k = 0; k<ctH->size(); k++) {
          CaloTowerRef t(ctH, k);
          if (&(**it) == &(*t)) {
            temp.push_back(k);
            break;
          }
        }
      }
      jet_calotwind->push_back(temp);
      
      if (jetTkAssColl.encode() != "") {
        std::vector<unsigned short> temp;       
        edm::Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
        iEvent.getByLabel(jetTkAssColl, jetTracksAssociation);
        
        for(reco::JetTracksAssociationCollection::const_iterator itass = jetTracksAssociation->begin(); itass != jetTracksAssociation->end(); ++itass) {
          if (&(*(itass->first)) != &(*j)) 
            continue;
          
          reco::TrackRefVector tracks = itass->second;
          
          jet_ntk[jet_n] = tracks.size();
          
          for (unsigned int ii = 0; ii < tracks.size(); ++ii) {
            for(unsigned int k = 0; k<tkH->size(); k++) {
              reco::TrackRef t(tkH, k);
              if (&(*(tracks[ii])) == &(*t) ) {
                temp.push_back(k);
                break;
              }
            }
          }
        }
        jet_tkind->push_back(temp);
      }
      
      jet_n++;
    }
  }
  
  
  if (strnome.find("PF",0) != std::string::npos){
    
    // take collections
    edm::Handle<reco::PFJetCollection> pfjetH;
    iEvent.getByLabel(jetColl, pfjetH);
    
    edm::Handle<CaloTowerCollection> ctH; 
    iEvent.getByLabel(calotowerColl, ctH);
    
    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(trackColl, tkH);
    
    edm::Handle<reco::VertexCollection> vtxH;
    iEvent.getByLabel(vertexColl, vtxH);
    
    jet_p4->Clear();
    jet_tkind->clear();
    jet_calotwind->clear();

    jet_n = 0;
    
    if (debug_level > 9)
      std::cout << "GlobeJets: Jet collection size: "<< pfjetH->size() << std::endl;
    
    if(iEvent.isRealData()){
      pfak5corr = pfak5corrdata;
    } else {
      pfak5corr = pfak5corrmc;
    }


    //PileupJetIdAlgo jetMVACalculator(*(jetMVAAlgos.begin()));
    PileupJetIdAlgo* jetMVACalculator = 0;
    if(algos_.size()>0) jetMVACalculator = algos_[0];
    // check if collection is present
    for(unsigned int i=0; i<pfjetH->size(); i++) {
      if (jet_n >= MAX_JETS) {
        std::cout << "GlobeJets: WARNING TOO MANY JETS: " << pfjetH->size() << " (allowed " << MAX_JETS << ")" << std::endl;
        break;
      }
      
      reco::PFJetRef j(pfjetH, i);
      
      // apply the cuts
      if (gCUT->cut(*j)) continue;
      // passed cuts
      
      reco::PFJet* correctedJet = (reco::PFJet*) j->clone();
      
      if (pfak5corr!="") {
	      const JetCorrector* corrector = JetCorrector::getJetCorrector(pfak5corr, iSetup);
        edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(pfjetH, i));
        jet_erescale[jet_n] = corrector->correction(*j, iEvent, iSetup);
      } else {
        jet_erescale[jet_n] = 1;
      }

      
      
      if (debug_level > 9 && jet_n<20) std::cout<<"pre "<<correctedJet->energy()<<std::endl;

      correctedJet->scaleEnergy(jet_erescale[jet_n]);

      if(correctedJet->pt() < 1) continue;
      
      pat::strbitset ret = pfLooseId.getBitTemplate();
      jet_pfloose[jet_n] = pfLooseId(*correctedJet, ret);
      
      //  scale the jets ref here

      new ((*jet_p4)[jet_n]) TLorentzVector();
      ((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(correctedJet->px(), correctedJet->py(), correctedJet->pz(), correctedJet->energy()); 
      jet_emfrac[jet_n] = correctedJet->chargedEmEnergyFraction() + correctedJet->neutralEmEnergyFraction() + correctedJet->chargedMuEnergyFraction();
      jet_hadfrac[jet_n] = correctedJet->chargedHadronEnergyFraction() + correctedJet->neutralHadronEnergyFraction();

      if (debug_level > 9 && jet_n<20) std::cout<<"post "<<correctedJet->energy()<<std::endl;

      if (debug_level > 9 && jet_n<20) std::cout<<"jet energy JECenergy "<<jet_n<<" "<<correctedJet->energy()<<" "<<jet_erescale[jet_n]<<std::endl;
      
      
      if (algos_.size()>0) {
        //compute id variables
        const reco::VertexCollection vertexCollection = *(vtxH.product()); 
        const reco::Vertex* selectedVtx  = &(*vertexCollection.begin());;
        const reco::Jet* thisjet = correctedJet;
        
        PileupJetIdentifier jetIdentifer_vars = jetMVACalculator->computeIdVariables( thisjet, jet_erescale[jet_n], selectedVtx, vertexCollection);
  
        jet_dRMean[jet_n]=jetIdentifer_vars.dRMean();
        jet_frac01[jet_n]=jetIdentifer_vars.frac01();
        jet_frac02[jet_n]=jetIdentifer_vars.frac02();
        jet_frac03[jet_n]=jetIdentifer_vars.frac03();
        jet_frac04[jet_n]=jetIdentifer_vars.frac04();
        jet_frac05[jet_n]=jetIdentifer_vars.frac05();
        jet_frac06[jet_n]=jetIdentifer_vars.frac06();
        jet_frac07[jet_n]=jetIdentifer_vars.frac07();
        jet_nNeutrals[jet_n]=jetIdentifer_vars.nNeutrals();
        jet_beta[jet_n]=jetIdentifer_vars.beta();
        jet_betaStar[jet_n]=jetIdentifer_vars.betaStar();
        jet_dZ[jet_n]=jetIdentifer_vars.dZ();
        jet_nCharged[jet_n]=jetIdentifer_vars.nCharged();
        jet_dR2Mean[jet_n]=jetIdentifer_vars.dR2Mean();
        jet_betaStarClassic[jet_n]=jetIdentifer_vars.betaStarClassic();
  
        for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
          PileupJetIdAlgo* ialgo = (algos_[imva]);
          ialgo->set(jetIdentifer_vars);
          PileupJetIdentifier id = ialgo->computeMva();
          mvas_[imva][jet_n] = id.mva() ;
          wp_levels_[imva][jet_n] = id.idFlag() ;
        }
      }

      float dy = 0;
      float dphi = 0;
      float dy_i = 0;
      float dphi_i = 0;
      float y_i = 0;
      float phi_i = 0;
      float dr_i = 0;
      float pt_i = 0;
      
      float pt_jt = j->pt();
      float y_jt = j->rapidity();
      float phi_jt = j->phi();
      unsigned int nCand = j->getPFConstituents().size();
      
      jet_pull_dy[jet_n]   = 0;
      jet_pull_dphi[jet_n] = 0;
      
      for(unsigned int c=0; c<nCand; c++) {
        reco::PFCandidatePtr iCand = j->getPFConstituent(c);
        y_i = iCand->rapidity(); //iCand rapidity
        dy_i = y_i - y_jt;
        
        phi_i = iCand->phi(); //iCand phi
        dphi_i = phi_i - phi_jt;
        
        dr_i=sqrt(dy_i*dy_i + phi_i*phi_i);
        
        pt_i= iCand->pt(); //iCand pt
        
        dy += pt_i*dr_i/pt_jt*dy_i;
        dphi += pt_i*dr_i/pt_jt*dphi_i;
      }
      
      jet_pull_dy[jet_n]   = dy;
      jet_pull_dphi[jet_n] = dphi;


          
      if (jetTkAssColl.encode() != "") {
        std::vector<unsigned short> temp;       
      
        edm::Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
        iEvent.getByLabel(jetTkAssColl, jetTracksAssociation);
        
        for(reco::JetTracksAssociationCollection::const_iterator itass = jetTracksAssociation->begin(); itass != jetTracksAssociation->end(); ++itass) {
          if (&(*(itass->first)) != &(*j)) 
            continue;
          
          reco::TrackRefVector tracks = itass->second;
          
          jet_ntk[jet_n] = tracks.size();
          
          for (unsigned int ii = 0; ii < tracks.size(); ++ii) {
      
            for(unsigned int k = 0; k<tkH->size(); k++) {
              reco::TrackRef t(tkH, k);
              if (&(*(tracks[ii])) == &(*t) ) {
                temp.push_back(k);
                break;
              }
            }
          }
        }
        jet_tkind->push_back(temp);
      }
      
      jet_n++;
      
    }
    
  }
  
  return true;
}
