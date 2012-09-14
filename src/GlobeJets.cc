#include <string>
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

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
  pfLooseId  = new PFJetIDSelectionFunctor( iConfig.getParameter<edm::ParameterSet>("pfLooseId") );
  std::string strnome = nome;
    
  mvas_.resize(jetMVAAlgos.size());
  wp_levels_.resize(jetMVAAlgos.size());
  mvas_ext_.resize(jetMVAAlgos.size());
  wp_levels_ext_.resize(jetMVAAlgos.size());
  algos_.resize(jetMVAAlgos.size());
  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
      algos_[imva] = new PileupJetIdAlgo((jetMVAAlgos.at(imva)));
  }


  if (strnome.find("PF",0) == std::string::npos){
    sprintf (a,"JetTrackAssociationColl_%s", nome);
    jetTkAssColl =  iConfig.getParameter<edm::InputTag>(a);
  }

  debug_level = iConfig.getParameter<int>("Debug_Level");
  jet_nvtx = iConfig.getParameter<unsigned int>("JetVertexToProcess");
  
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
  
  sprintf(a1, "jet_%s_area", nome);
  sprintf(a2, "jet_%s_area[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_area, a2);

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
  sprintf(a2, "jet_%s_nNeutrals[jet_%s_n]/I", nome, nome);
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
  sprintf(a2, "jet_%s_nCharged[jet_%s_n]/I", nome, nome);
  tree->Branch(a1, &jet_nCharged, a2);

  sprintf(a1, "jet_%s_rmsCand", nome);
  sprintf(a2, "jet_%s_rmsCand[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_rmsCand, a2);

  sprintf(a1, "jet_%s_ptD", nome);
  sprintf(a2, "jet_%s_ptD[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_ptD, a2);

  sprintf(a1, "jet_%s_axis1", nome);
  sprintf(a2, "jet_%s_axis1[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_axis1, a2);

  sprintf(a1, "jet_%s_axis2", nome);
  sprintf(a2, "jet_%s_axis2[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_axis2, a2);

  sprintf(a1, "jet_%s_pull", nome);
  sprintf(a2, "jet_%s_pull[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_pull, a2);

  sprintf(a1, "jet_%s_tana", nome);
  sprintf(a2, "jet_%s_tana[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_tana, a2);

  sprintf(a1, "jet_%s_rmsCand_QC", nome);
  sprintf(a2, "jet_%s_rmsCand_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_rmsCand_QC, a2);

  sprintf(a1, "jet_%s_ptD_QC", nome);
  sprintf(a2, "jet_%s_ptD_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_ptD_QC, a2);

  sprintf(a1, "jet_%s_axis1_QC", nome);
  sprintf(a2, "jet_%s_axis1_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_axis1_QC, a2);

  sprintf(a1, "jet_%s_axis2_QC", nome);
  sprintf(a2, "jet_%s_axis2_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_axis2_QC, a2);

  sprintf(a1, "jet_%s_pull_QC", nome);
  sprintf(a2, "jet_%s_pull_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_pull_QC, a2);

  sprintf(a1, "jet_%s_tana_QC", nome);
  sprintf(a2, "jet_%s_tana_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_tana_QC, a2);

  sprintf(a1, "jet_%s_Rchg", nome);
  sprintf(a2, "jet_%s_Rchg[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_Rchg, a2);

  sprintf(a1, "jet_%s_Rneutral", nome);
  sprintf(a2, "jet_%s_Rneutral[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_Rneutral, a2);

  sprintf(a1, "jet_%s_R", nome);
  sprintf(a2, "jet_%s_R[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_R, a2);

  sprintf(a1, "jet_%s_Rchg_QC", nome);
  sprintf(a2, "jet_%s_Rchg_QC[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_Rchg_QC, a2);


  sprintf(a1, "jet_%s_dR2Mean", nome);
  sprintf(a2, "jet_%s_dR2Mean[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_dR2Mean, a2);

  sprintf(a1, "jet_%s_betaStarClassic", nome);
  sprintf(a2, "jet_%s_betaStarClassic[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_betaStarClassic, a2);

  sprintf(a1, "jet_%s_beta_ext", nome);
  tree->Branch(a1, "std::vector<std::vector<float> >", &jet_beta_ext);
  sprintf(a1, "jet_%s_betaStar_ext", nome);
  tree->Branch(a1, "std::vector<std::vector<float> >", &jet_betaStar_ext);
  sprintf(a1, "jet_%s_betaStarClassic_ext", nome);
  tree->Branch(a1, "std::vector<std::vector<float> >", &jet_betaStarClassic_ext);

  sprintf(a1, "jet_%s_nvtx", nome);
  sprintf(a2, "jet_%s_nvtx/i", nome);
  tree->Branch(a1, &jet_nvtx, a2);
  
  sprintf(a1, "jet_%s_pfloose", nome);
  sprintf(a2, "jet_%s_pfloose[jet_%s_n]/O", nome, nome);
  tree->Branch(a1, &jet_pfloose, a2);

  sprintf(a1, "jet_%s_csvBtag", nome);
  sprintf(a2, "jet_%s_csvBtag[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_csvBtag, a2);

  sprintf(a1, "jet_%s_csvMvaBtag", nome);
  sprintf(a2, "jet_%s_csvMvaBtag[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_csvMvaBtag, a2);

  sprintf(a1, "jet_%s_jetProbBtag", nome);
  sprintf(a2, "jet_%s_jetProbBtag[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_jetProbBtag, a2);

  sprintf(a1, "jet_%s_tcheBtag", nome);
  sprintf(a2, "jet_%s_tcheBtag[jet_%s_n]/F", nome, nome);
  tree->Branch(a1, &jet_tcheBtag, a2);

 
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

      mvas_ext_[imva] = new std::vector<std::vector<float> >();
      wp_levels_ext_[imva] = new std::vector<std::vector<int> >();

      sprintf(a1, "jet_%s_%s_mva_ext", nome, mvalabel.c_str());
      tree->Branch(a1, "std::vector<std::vector<float> >", &mvas_ext_[imva]);

      sprintf(a1, "jet_%s_%s_wp_level_ext", nome, mvalabel.c_str());
      tree->Branch(a1, "std::vector<std::vector<int> >", &wp_levels_ext_[imva]);
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
    
    jet_p4->Delete();
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
      //if (gCUT->cut(*j)) continue;
      // passed cuts
      
      new ((*jet_p4)[jet_n]) TLorentzVector();
      ((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(j->px(), j->py(), j->pz(), j->energy()); 
      jet_area[jet_n] = j->jetArea();
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
    
    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexBJetTags;
    if( strnome=="algoPF1" )
      iEvent.getByLabel("newPFCombinedSecondaryVertexBPFJetTags", combinedSecondaryVertexBJetTags);
    else
      iEvent.getByLabel("newPFchsCombinedSecondaryVertexBPFJetTags", combinedSecondaryVertexBJetTags);
    
    edm::Handle<reco::JetTagCollection> combinedSecondaryVertexMVABJetTags;
    if( strnome=="algoPF1" )
      iEvent.getByLabel("newPFCombinedSecondaryVertexMVABPFJetTags", combinedSecondaryVertexMVABJetTags);
    else
      iEvent.getByLabel("newPFchsCombinedSecondaryVertexMVABPFJetTags", combinedSecondaryVertexMVABJetTags);
    
    edm::Handle<reco::JetTagCollection> jetProbabilityBJetTags;
    if( strnome=="algoPF1" )
      iEvent.getByLabel("newPFJetProbabilityBPFJetTags", jetProbabilityBJetTags);
    else
      iEvent.getByLabel("newPFchsJetProbabilityBPFJetTags", jetProbabilityBJetTags);
    
    edm::Handle<reco::JetTagCollection> trackCountingHighEffBJetTags;
    if( strnome=="algoPF1" )
      iEvent.getByLabel("newPFTrackCountingHighEffBJetTags", trackCountingHighEffBJetTags);
    else
      iEvent.getByLabel("newPFchsTrackCountingHighEffBJetTags", trackCountingHighEffBJetTags);

    
    jet_p4->Delete();
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


    PileupJetIdAlgo* jetMVACalculator = 0;
    if(algos_.size()>0) jetMVACalculator = algos_[0];
    size_t n_jets = std::min(pfjetH->size(),(size_t)MAX_JETS);
    size_t n_vtx = std::min(vtxH->size(), (size_t)jet_nvtx); //vtxH->size();
    for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
      mvas_ext_[imva]->clear();
      mvas_ext_[imva]->resize(n_jets, std::vector<float>(n_vtx,-999.));
      wp_levels_ext_[imva]->clear();
      wp_levels_ext_[imva]->resize(n_jets, std::vector<int>(n_vtx,-999.));
      jet_beta_ext.clear();
      jet_beta_ext.resize(n_jets,std::vector<float>(n_vtx,-999.));
      jet_betaStar_ext.clear();
      jet_betaStar_ext.resize(n_jets,std::vector<float>(n_vtx,-999.));
      jet_betaStarClassic_ext.clear();
      jet_betaStarClassic_ext.resize(n_jets,std::vector<float>(n_vtx,-999.));
      
    }

    // check if collection is present
    for(unsigned int i=0; i<pfjetH->size(); i++) {
      if (jet_n >= MAX_JETS) {
	std::cout << "GlobeJets: WARNING TOO MANY JETS: " << pfjetH->size() << " (allowed " << MAX_JETS << ")" << std::endl;
	break;
      }
      
      reco::PFJetRef j(pfjetH, i);
      
      // apply the cuts
      if (gCUT->cut(*j)) 
	continue;
      // passed cuts
      
      reco::PFJet* correctedJet = (reco::PFJet*) j->clone();
      
      if (pfak5corr!="") {
	const JetCorrector* corrector = JetCorrector::getJetCorrector(pfak5corr, iSetup);
        edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::PFJetCollection>(pfjetH, i));
        jet_erescale[jet_n] = corrector->correction(*j, iEvent, iSetup);
	//delete corrector;
	
      } else {
        jet_erescale[jet_n] = 1;
      }
      
      if (debug_level > 9 && jet_n<20) std::cout<<"pre "<<correctedJet->energy()<<std::endl;

      correctedJet->scaleEnergy(jet_erescale[jet_n]);

      //if(correctedJet->pt() < 1) continue;
      
      pat::strbitset ret = (*pfLooseId).getBitTemplate();
      jet_pfloose[jet_n] = (*pfLooseId)(*j, ret);
      
      //  scale the jets ref here

      new ((*jet_p4)[jet_n]) TLorentzVector();
      ((TLorentzVector *)jet_p4->At(jet_n))->SetXYZT(correctedJet->px(), correctedJet->py(), correctedJet->pz(), correctedJet->energy()); 
      jet_area[jet_n] = correctedJet->jetArea();
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

	size_t n_vtx = std::min(vtxH->size(), (size_t)jet_nvtx); //vtxH->size();
	//for(size_t vtx=0; vtx<vertexCollection.size(); ++vtx) {
	for(size_t vtx=0; vtx<n_vtx; ++vtx) {
	  PileupJetIdentifier ext_vars = jetMVACalculator->computeIdVariables( thisjet, jet_erescale[jet_n], &vertexCollection[vtx], vertexCollection);
	  jet_beta_ext[jet_n][vtx]=ext_vars.beta();
	  jet_betaStar_ext[jet_n][vtx]=ext_vars.betaStar();
	  jet_betaStarClassic_ext[jet_n][vtx]=ext_vars.betaStarClassic();
	  
	  for(unsigned int imva=0; imva<jetMVAAlgos.size(); imva++){
	    PileupJetIdAlgo* ialgo = (algos_[imva]);
	    ialgo->set(ext_vars);
	    PileupJetIdentifier id = ialgo->computeMva();
	    (*mvas_ext_[imva])[jet_n][vtx] = id.mva() ;
	    (*wp_levels_ext_[imva])[jet_n][vtx] = id.idFlag();
	  }
	}
      }

      delete correctedJet;

      float dy = 0;
      float dphi = 0;
      float dy_i = 0;
      float dphi_i = 0;
      float deta_i = 0;
      float y_i = 0;
      float phi_i = 0;
      float eta_i = 0;
      float dr_i = 0;
      float pt_i = 0;
      
      float pt_jt = j->pt();
      float y_jt = j->rapidity();
      float phi_jt = j->phi();
      float eta_jt = j->eta();
      unsigned int nCand = j->getPFConstituents().size();
      
      
      // Quark-Gluon discrimination variable computation --- BEGIN

      jet_pull_dy[jet_n]   = 0;
      jet_pull_dphi[jet_n] = 0;

      jet_ptD[jet_n] = -1.;
      jet_rmsCand[jet_n] = -1.;
      jet_axis1[jet_n] = -1.;
      jet_axis2[jet_n] = -1.;
      jet_pull[jet_n] = -1.;
      jet_tana[jet_n] = -1.;

      jet_ptD_QC[jet_n] = -1.;
      jet_rmsCand_QC[jet_n] = -1.;
      jet_axis1_QC[jet_n] = -1.;
      jet_axis2_QC[jet_n] = -1.;
      jet_pull_QC[jet_n] = -1.;
      jet_tana_QC[jet_n] = -1.;

      jet_Rchg[jet_n] = -1.;
      jet_Rneutral[jet_n] = -1.;
      jet_R[jet_n] = -1.;
      jet_Rchg_QC[jet_n] = -1.;

      //float rmsCand=  -999.;
      //float ptD=      -999.;
      //float axis1=    -999.;
      //float axis2=    -999.;
      //float pull=     -999.;
      //float tana =    -999.;
      //float rmsCand_QC=  -999.;
      //float ptD_QC=      -999.;
      //float axis1_QC=    -999.;
      //float axis2_QC=    -999.;
      //float pull_QC=     -999.;
      //float tana_QC =    -999.;
      
      float SumW=0;
      float SumW2=0;
      float SumDeta=0;
      float SumDeta2=0;
      float SumDphi=0;
      float SumDphi2=0;
      float SumDetaDphi=0;
      
      float SumW_QC=0;
      float SumW2_QC=0;
      float SumDeta_QC=0;
      float SumDeta2_QC=0;
      float SumDphi_QC=0;
      float SumDphi2_QC=0;
      float SumDetaDphi_QC=0;

      float pTMax(0.0),pTMaxChg(0.0),pTMaxNeutral(0.0),pTMaxChg_QC(0.0);
      
      std::vector<bool> jetPart_forMult,jetPart_forAxis;

      
      
      // first loop on cands

      for(unsigned int c=0; c<nCand; c++) {


        reco::PFCandidatePtr iCand = j->getPFConstituent(c);
        y_i = iCand->rapidity(); //iCand rapidity
        dy_i = y_i - y_jt;

        eta_i = iCand->eta(); //iCand rapidity
        deta_i = eta_i - eta_jt;
        
        phi_i = iCand->phi(); //iCand phi
        dphi_i = phi_i - phi_jt;
        double dphi_i_wrap = 2*atan(tan((phi_i-phi_jt)/2));      
        
        dr_i=sqrt(dy_i*dy_i + phi_i*phi_i);
        
        pt_i= iCand->pt(); //iCand pt
        
        dy += pt_i*dr_i/pt_jt*dy_i;
        dphi += pt_i*dr_i/pt_jt*dphi_i;


        reco::TrackRef itrk ;
        if (iCand.isNonnull())
          itrk = (*iCand).trackRef();
        if (pt_i > pTMax) 
          pTMax = pt_i;
        if (itrk.isNonnull() && pt_i > pTMaxChg) 
          pTMaxChg = pt_i;
        if (!itrk.isNonnull() && pt_i > pTMaxNeutral) 
          pTMaxNeutral = pt_i;

        
        bool trkForAxis = false;
        bool trkForMult = false;
        
        //-----matching with vertex tracks-------
        if (!itrk.isNonnull()) { 
          trkForMult = true;
          trkForAxis = true;
        }
        else {
          float dZmin = 999;
          int index_min = 999;
          reco::VertexCollection::const_iterator vtxClose;
          for(unsigned ivtx = 0;ivtx < n_vtx;ivtx++) {
            float dZ_cut = fabs(itrk->dz((*vtxH)[ivtx].position()));
            float sumpT = 0;
            for(reco::Vertex::trackRef_iterator itk = (*vtxH)[ivtx].tracks_begin();itk!=(*vtxH)[ivtx].tracks_end(); ++itk) {
              sumpT = sumpT + ((*itk)->pt())*((*itk)->pt());
            }
            if (dZ_cut < dZmin) {
              dZmin = dZ_cut;
              index_min = ivtx;
            }
          }//Loop over vertices 

          if (index_min == 0) {
            float dz = itrk->dz((*vtxH)[0].position());
            float d0 = itrk->dxy((*vtxH)[0].position());
            float vtx_xError = (*vtxH)[0].xError();
            float vtx_yError = (*vtxH)[0].yError();
            float vtx_zError = (*vtxH)[0].zError();
            float d0_sigma=sqrt(pow(itrk->d0Error(),2) + pow(vtx_xError,2) + pow(vtx_yError,2));
            float dz_sigma=sqrt(pow(itrk->dzError(),2) + pow(vtx_zError,2));
            if (itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.) {
              trkForAxis = true;
              if (fabs(d0/d0_sigma) < 5.)
                trkForMult = true;
            }
          }
          if (pt_i > pTMaxChg_QC && trkForAxis)  
            pTMaxChg_QC = pt_i;
        }// for charged particles only


        jetPart_forMult.push_back(trkForMult);
        jetPart_forAxis.push_back(trkForAxis);
    
        SumW+=pt_i;
        SumW2+=pt_i*pt_i;
        SumDeta+=pt_i*pt_i*deta_i;
        SumDeta2+=pt_i*pt_i*deta_i*deta_i;
        SumDphi+=pt_i*pt_i*dphi_i_wrap;
        SumDphi2+=pt_i*pt_i*dphi_i_wrap*dphi_i_wrap;
        SumDetaDphi+=pt_i*pt_i*deta_i*dphi_i_wrap;
        if (trkForAxis) {
          SumW_QC+=pt_i;
          SumW2_QC+=pt_i*pt_i;
          SumDeta_QC+=pt_i*pt_i*deta_i;
          SumDeta2_QC+=pt_i*pt_i*deta_i*deta_i;
          SumDphi_QC+=pt_i*pt_i*dphi_i_wrap;
          SumDphi2_QC+=pt_i*pt_i*dphi_i_wrap*dphi_i_wrap;
          SumDetaDphi_QC+=pt_i*pt_i*deta_i*dphi_i_wrap;
        }

      } //first loop on PFCandidates


      float ave_deta = SumDeta/SumW2;
      float ave_dphi = SumDphi/SumW2;
      float ave_deta2 = SumDeta2/SumW2;
      float ave_dphi2 = SumDphi2/SumW2;
      float a = ave_deta2-ave_deta*ave_deta;
      float b = ave_dphi2-ave_dphi*ave_dphi;
      float c = -(SumDetaDphi/SumW2-ave_deta*ave_dphi);
      float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
      if (a+b+delta >= 0) {
        jet_axis1[jet_n] = sqrt(0.5*(a+b+delta));
      }
      if (a+b-delta >= 0) {  
        jet_axis2[jet_n] = sqrt(0.5*(a+b-delta));
      }
      if (c != 0) {
        jet_tana[jet_n] = 0.5*(b-a+delta)/c;
      }

      jet_ptD[jet_n] = sqrt( SumW2/ (SumW*SumW));


      float ave_deta_QC = SumDeta_QC/SumW2_QC;
      float ave_dphi_QC = SumDphi_QC/SumW2_QC;
      float ave_deta2_QC = SumDeta2_QC/SumW2_QC;
      float ave_dphi2_QC = SumDphi2_QC/SumW2_QC;
      float a_QC = ave_deta2_QC-ave_deta_QC*ave_deta_QC;
      float b_QC = ave_dphi2_QC-ave_dphi_QC*ave_dphi_QC;
      float c_QC = -(SumDetaDphi_QC/SumW2_QC-ave_deta_QC*ave_dphi_QC);
      float delta_QC = sqrt(fabs((a_QC-b_QC)*(a_QC-b_QC)+4*c_QC*c_QC));
      if (a_QC+b_QC+delta_QC >= 0) {
        jet_axis1_QC[jet_n] = sqrt(0.5*(a_QC+b_QC+delta_QC));
      }
      if (a_QC+b_QC-delta_QC >= 0) {  
        jet_axis2_QC[jet_n] = sqrt(0.5*(a_QC+b_QC-delta_QC));
      }
      if (c != 0) {
        jet_tana_QC[jet_n] = 0.5*(b-a+delta)/c;
      }

      jet_ptD_QC[jet_n] = sqrt( SumW2_QC/ (SumW_QC*SumW_QC));



      //-------second loop to calculate higher moments

      float ddetaR_sum(0.0), ddphiR_sum(0.0),ddetaR_sum_QC(0.0), ddphiR_sum_QC(0.0);
      float sum_ddR = 0.;
      float sum_ddR_QC = 0.;

      for(int icand=0; icand<j->nConstituents(); ++icand) {

        double pt=j->getJetConstituentsQuick()[icand]->pt();
        double eta=j->getJetConstituentsQuick()[icand]->eta();
        double phi=j->getJetConstituentsQuick()[icand]->phi();
        double dphi = 2*atan(tan((phi-phi_jt)/2));      
        double deta = eta-eta_jt;
        float weight = pt*pt;
        float ddeta, ddphi,ddR;
        ddeta = deta - ave_deta ;
        ddphi = 2*atan(tan(( dphi - ave_dphi)/2.)) ;
        ddR = sqrt(ddeta*ddeta + ddphi*ddphi);
        sum_ddR += ddR *ddR* weight;
        ddetaR_sum += ddR*ddeta*weight;
        ddphiR_sum += ddR*ddphi*weight;
        if (jetPart_forAxis[icand]) { 
          float ddeta_QC = deta - ave_deta_QC ;
          float ddphi_QC = 2*atan(tan(( dphi - ave_dphi_QC)/2.)) ;
          float ddR_QC = sqrt(ddeta_QC*ddeta_QC + ddphi_QC*ddphi_QC);
          sum_ddR_QC += ddR_QC *ddR_QC* weight;
          ddetaR_sum_QC += ddR_QC*ddeta_QC*weight;
          ddphiR_sum_QC += ddR_QC*ddphi_QC*weight;
        }
     
      }//second loop over constituents  



      if (SumW2 > 0) {
        float ddetaR_ave = ddetaR_sum/SumW2;
        float ddphiR_ave = ddphiR_sum/SumW2;
        jet_pull[jet_n] = sqrt(ddetaR_ave*ddetaR_ave+ddphiR_ave*ddphiR_ave);
      }

      if (SumW2_QC > 0) {
        float ddetaR_ave_QC = ddetaR_sum_QC/SumW2_QC;
        float ddphiR_ave_QC = ddphiR_sum_QC/SumW2_QC;
        jet_pull_QC[jet_n] = sqrt(ddetaR_ave_QC*ddetaR_ave_QC+ddphiR_ave_QC*ddphiR_ave_QC);
      }


      jet_rmsCand[jet_n] = sqrt( sum_ddR / SumW2);
      jet_rmsCand_QC[jet_n] = sqrt( sum_ddR_QC / SumW2_QC);

      jet_Rchg[jet_n] = pTMaxChg/SumW;
      jet_Rneutral[jet_n] = pTMaxNeutral/SumW;
      jet_R[jet_n] = pTMax/SumW;
      jet_Rchg_QC[jet_n] = pTMaxChg_QC/SumW_QC;


      
      jet_pull_dy[jet_n]   = dy;
      jet_pull_dphi[jet_n] = dphi;


      // Quark-Gluon discrimination variable computation --- END



      std::vector<unsigned short> temp;       
      reco::TrackRefVector tracks = j->getTrackRefs();
      jet_ntk[jet_n] = tracks.size();
          
      for (unsigned int ii = 0; ii < tracks.size(); ++ii) {
        for(unsigned int k = 0; k<tkH->size(); k++) {
          reco::TrackRef t(tkH, k);
          if(gCUT->cut(*t)) continue; 
          if (&(*(tracks[ii])) == &(*t) ) {
            temp.push_back(k);
            break;
          }
        }
      }
      jet_tkind->push_back(temp);


      // btags:
      jet_csvBtag    [jet_n] =  (*combinedSecondaryVertexBJetTags)[i].second ;
      jet_csvMvaBtag [jet_n] =  (*combinedSecondaryVertexMVABJetTags)[i].second ;
      jet_jetProbBtag[jet_n] =  (*jetProbabilityBJetTags)[i].second ;
      jet_tcheBtag   [jet_n] =  (*trackCountingHighEffBJetTags)[i].second ;

      
      jet_n++;
      
    }
    
  }
  
  return true;
}
