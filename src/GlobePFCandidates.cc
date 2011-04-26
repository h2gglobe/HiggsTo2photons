#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePFCandidates.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <iostream>

GlobePFCandidates::GlobePFCandidates(const edm::ParameterSet& iConfig) {
  
  pfColl = iConfig.getParameter<edm::InputTag>("PFCandidateColl");
  photonCollStd =  iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  PFIsoOuterConeSize = iConfig.getParameter<double>("PFIsoOuterCone");
  debug_level = iConfig.getParameter<int>("Debug_Level");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobePFCandidates::defineBranch(TTree* tree) {
  
  pfcand_p4 = new TClonesArray("TLorentzVector", MAX_PFCANDS);
  pfcand_poscalo = new TClonesArray("TVector3", MAX_PFCANDS);
  
  tree->Branch("pfcand_p4", "TClonesArray", &pfcand_p4, 32000, 0);
  tree->Branch("pfcand_poscalo", "TClonesArray", &pfcand_poscalo, 32000, 0);

  tree->Branch("pfcand_n", &pfcand_n, "pfcand_n/I");
  tree->Branch("pfcand_pdgid",&pfcand_pdgid,"pfcand_pdgid[pfcand_n]/I");
  tree->Branch("pfcand_tkind", &pfcand_tkind, "pfcand_tkind[pfcand_n]/I");
  //tree->Branch("pfcand_gsfind", &pfcand_gsfind, "pfcand_gsfind[pfcand_n]/I");
  //tree->Branch("pfcand_muind", &pfcand_muind, "pfcand_muind[pfcand_n]/I");
  tree->Branch("pfcand_ecalenergy", &pfcand_ecalEnergy, "pfcand_ecalenergy[pfcand_n]/F");
  tree->Branch("pfcand_hcalenergy", &pfcand_hcalEnergy, "pfcand_hcalenergy[pfcand_n]/F");
  tree->Branch("pfcand_rawecalenergy", &pfcand_rawEcalEnergy, "pfcand_rawecalenergy[pfcand_n]/F");
  tree->Branch("pfcand_rawhcalenergy", &pfcand_rawHcalEnergy, "pfcand_rawhcalenergy[pfcand_n]/F");
  tree->Branch("pfcand_ps1energy", &pfcand_ps1Energy, "pfcand_ps1energy[pfcand_n]/F");
  tree->Branch("pfcand_ps2energy", &pfcand_ps2Energy, "pfcand_ps2energy[pfcand_n]/F");
  tree->Branch("pfcand_momerr", &pfcand_momErr, "pfcand_momerr[pfcand_n]/F");
  tree->Branch("pfcand_mva_e_pi", &pfcand_mva_e_pi, "pfcand_mva_e_pi[pfcand_n]/F");
  tree->Branch("pfcand_mva_e_mu", &pfcand_mva_e_mu, "pfcand_mva_e_mu[pfcand_n]/F");
  tree->Branch("pfcand_mva_pi_mu", &pfcand_mva_pi_mu, "pfcand_mva_pi_mu[pfcand_n]/F");
  tree->Branch("pfcand_mva_nothing_gamma", &pfcand_mva_nothing_gamma, "pfcand_mva_nothing_gamma[pfcand_n]/F");
  tree->Branch("pfcand_mva_nothing_nh", &pfcand_mva_nothing_nh, "pfcand_mva_nothing_nh[pfcand_n]/F");
  tree->Branch("pfcand_mva_gamma_nh", &pfcand_mva_gamma_nh, "pfcand_mva_gamma_nh[pfcand_n]/F");
  tree->Branch("pfcand_vz",&pfcand_vz,"pfcand_vz[pfcand_n]/F");
  tree->Branch("pfcand_overlappho",&pfcand_overlappho,"pfcand_overlappho[pfcand_n]/i");
  //tree->Branch("pfcand_overlappho",&pfcand_overlappho,"pfcand_overlappho[pfcand_n][pho_n]/I");
}

bool GlobePFCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeTracks* tks, GlobeMuons* mus, GlobePhotons* phos) {

  // take collections
  edm::Handle<reco::PFCandidateCollection> pfH;
  iEvent.getByLabel(pfColl, pfH);

  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);


  pfcand_p4->Clear(); 
  pfcand_poscalo->Clear(); 
  pfcand_n = 0;
  //pfcandtimespho_n = 0;

  if (debug_level > 9)
    std::cout << "GlobePFCandidates: PFCandidates collection size: "<< pfH->size() << std::endl;
  
  for(reco::PFCandidateCollection::const_iterator ipf = pfH->begin(); ipf != pfH->end(); ipf++) {

    if (pfcand_n >= MAX_PFCANDS) {
      std::cout << "GlobePFCandidates: WARNING TOO MANY CANDIDATES: " << pfH->size() << " (allowed " << MAX_PFCANDS << ")" << std::endl;
      break;
    }
    
    //reco::PFCandidate pfc = reco::PFCandidate(ipf);
    
    if(gCUT->cut(*ipf))
      continue;

    bool keepPfCand = false;
      
    if (phos!=0){
      for (Int_t i=0; i<phos->pho_n; i++){
	TVector3* PhoEcalPos = (TVector3*)phos->pho_calopos->At(i);
	TVector3* PfCandEcalPos = new TVector3(ipf->positionAtECALEntrance().x(), 
					       ipf->positionAtECALEntrance().y(), 
					       ipf->positionAtECALEntrance().z());
	if (PfCandEcalPos->X()!=0 && PfCandEcalPos->Y()!=0 && PfCandEcalPos->Z()!=0){
	  double dEta = PhoEcalPos->Eta() - PfCandEcalPos->Eta();
	  double dPhi = PhoEcalPos->Phi() - PfCandEcalPos->Phi();
	  
	  if (sqrt(dEta*dEta+dPhi*dPhi)<PFIsoOuterConeSize) keepPfCand=true;
	}
      }
    }
    if (keepPfCand==false) continue;


    //pfcand_pdgid[pfcand_n] = ipf->translateTypeToPdgId(ipf->particleId());
    pfcand_pdgid[pfcand_n] = ipf->pdgId();
    pfcand_ecalEnergy[pfcand_n] = ipf->ecalEnergy();
    pfcand_hcalEnergy[pfcand_n] = ipf->hcalEnergy();
    pfcand_rawEcalEnergy[pfcand_n] = ipf->rawEcalEnergy();
    pfcand_rawHcalEnergy[pfcand_n] = ipf->rawHcalEnergy();
    pfcand_ps1Energy[pfcand_n] = ipf->pS1Energy();
    pfcand_ps2Energy[pfcand_n] = ipf->pS2Energy();
    pfcand_momErr[pfcand_n] = ipf->deltaP();
    pfcand_mva_e_pi[pfcand_n] = ipf->mva_e_pi();
    pfcand_mva_e_mu[pfcand_n] = ipf->mva_e_mu();
    pfcand_mva_pi_mu[pfcand_n] = ipf->mva_pi_mu();
    pfcand_mva_nothing_gamma[pfcand_n] = ipf->mva_nothing_gamma();
    pfcand_mva_nothing_nh[pfcand_n] = ipf->mva_nothing_nh();
    pfcand_mva_gamma_nh[pfcand_n] = ipf->mva_gamma_nh();
    pfcand_vz[pfcand_n] = ipf->vz();

    new ((*pfcand_p4)[pfcand_n]) TLorentzVector();
    ((TLorentzVector *)pfcand_p4->At(pfcand_n))->SetXYZT(ipf->px(), ipf->py(), ipf->pz(), ipf->energy());
     
    new ((*pfcand_poscalo)[pfcand_n]) TVector3();
    ((TVector3 *)pfcand_poscalo->At(pfcand_n))->SetXYZ(ipf->positionAtECALEntrance().x(), 
						       ipf->positionAtECALEntrance().y(), 
						       ipf->positionAtECALEntrance().z());

    pfcand_tkind[pfcand_n] = -1;  
    if (tks != 0 && ipf->trackRef().isNonnull()) {
      TVector3 vTk = TVector3(ipf->trackRef()->px(), 
                              ipf->trackRef()->py(), 
                              ipf->trackRef()->pz());
      
      float drmin = 0.2;
      for (Int_t i=0; i<tks->tk_n; i++) {
        TVector3 temp = ((TLorentzVector*)tks->tk_p4->At(i))->Vect();
        float dr = temp.DeltaR(vTk);
        if (dr < drmin) {
          drmin = dr;
          pfcand_tkind[pfcand_n] = i;
        }
      }
    }

    pfcand_muind[pfcand_n] = -1;
    if (mus != 0 && ipf->muonRef().isNonnull()) {
      TVector3 vMu = TVector3(ipf->muonRef()->px(), 
                              ipf->muonRef()->py(), 
                              ipf->muonRef()->pz());

      float drmin = 0.2;
      for (Int_t i=0; i<mus->mu_n; i++) {
        TLorentzVector* temp = (TLorentzVector*)mus->mu_p4->At(i);
        float dr = temp->Vect().DeltaR(vMu);
        if (dr < drmin) {
          drmin = dr;
          pfcand_muind[pfcand_n] = i;
        }
      }
    }
      

    
    if (phos!=0){

      int coding = 1;
      
      //pho_n = 0;
      pfcand_overlappho[pfcand_n] = 0;

      /*
      //pfcand_overlappho[pfcand_n][pho_n] = 0;

      for (Int_t i=0; i<phos->pho_n; i++){
	TLorentzVector* temp = (TLorentzVector*)phos->pho_p4->At(i);
	for( reco::PhotonCollection::const_iterator  iPho = phoH->begin(); iPho != phoH->end(); iPho++) {

	  if (temp->Px()==iPho->px() && temp->Py()==iPho->py() && temp->Pz()==iPho->pz() && temp->Energy()==iPho->energy()){
	    
	    cout << "Photon egal"<<endl;

	    //get supercluster
	    const reco::SuperCluster* supercluster = new reco::SuperCluster(*(iPho->superCluster()));
	    std::vector<const reco::SuperCluster*> sc;
	    sc.push_back(supercluster);

	    //get pfcandidate clusters
	    for(unsigned iele=0; iele<ipf->elementsInBlocks().size(); ++iele) {

	      reco::PFBlockRef blockRef = ipf->elementsInBlocks()[iele].first;
	      unsigned elementIndex = ipf->elementsInBlocks()[iele].second;

	      if(!blockRef.isNull()){

		const edm::OwnVector< reco::PFBlockElement >&  elements = (*blockRef).elements();
		const reco::PFBlockElement & pfbe (elements[elementIndex]);

		if(pfbe.type()==reco::PFBlockElement::ECAL || pfbe.type()==reco::PFBlockElement::PS1 || pfbe.type()==reco::PFBlockElement::PS2 || pfbe.type()==reco::PFBlockElement::HCAL || pfbe.type()==reco::PFBlockElement::SC || pfbe.type()==reco::PFBlockElement::BREM){

		  reco::PFClusterRef myPFClusterRef = pfbe.clusterRef();
		  if(!myPFClusterRef.isNull()){
		    
		    const reco::PFCluster & myPFCluster (*myPFClusterRef);
		    //cout << "PFcand has a ClusterRef"<<endl;
		    int hasOverlap = ClusterClusterMapping::checkOverlap(myPFCluster,sc);
		    //if (hasOverlap==-1) cout << "NO Overlap with E/gamma SC"<<endl;
		    //else cout << "Overlap with E/gamma SC"<<endl;
		    
		    //if (hasOverlap!=-1) pfcand_overlappho[pfcandtimespho_n] = 1;
		    if (hasOverlap!=-1) {
		      //pfcand_overlappho[pfcand_n][pho_n] = 1;
		      pfcand_overlappho[pfcand_n] = (pfcand_overlappho[pfcand_n] | coding);
		    }
		    //cout << "pfcand_overlappho[pfcand_n][pho_n]=" <<pfcand_overlappho[pfcand_n][pho_n]<<endl;

		  }
		}
	      }
	    }
	    //pho_n++;
	    coding = coding*2;

	  }

	}
      }
      */
    }
    
    

    pfcand_n++;
  }
     
  return true;
}

