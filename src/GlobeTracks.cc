#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"

//----------------------------------------------------------------------

GlobeTracks::GlobeTracks(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  doAodSim     = iConfig.getParameter<bool>("doAodSim");
  trackColl    = iConfig.getParameter<edm::InputTag>("TrackColl");
  trackColl2   = iConfig.getParameter<edm::InputTag>("TrackColl2");
  tpColl       = iConfig.getParameter<edm::InputTag>("TPColl");
  assocLabel   = iConfig.getParameter<std::string>("AssocLabel");
  debug_level  = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

//----------------------------------------------------------------------

void GlobeTracks::defineBranch(TTree* tree) {
  
  tk_p4 = new TClonesArray("TLorentzVector", MAX_TRACKS);
  tk_vtx_pos = new TClonesArray("TVector3", MAX_TRACKS);
  
  tree->Branch("tk_n", &tk_n, "tk_n/I");
  tree->Branch("tk_p4", "TClonesArray", &tk_p4, 32000, 0);
  tree->Branch("tk_vtx_pos", "TClonesArray", &tk_vtx_pos, 32000, 0);
  tree->Branch("tk_nhits", &tk_nhits, "tk_nhits[tk_n]/I");
  tree->Branch("tk_charge", &tk_charge, "tk_charge[tk_n]/I");
  tree->Branch("tk_nlosthit", &tk_nlosthit,"tk_nlosthit[tk_n]/I" );
  tree->Branch("tk_tpind", &tk_tpind,"tk_tpind[tk_n]/I" );
  tree->Branch("tk_chi2", &tk_chi2,"tk_chi2[tk_n]/F" );
  tree->Branch("tk_dof", &tk_dof, "tk_dof[tk_n]/F");
  tree->Branch("tk_d0", &tk_d0, "tk_d0[tk_n]/F");
  tree->Branch("tk_dz", &tk_dz, "tk_dz[tk_n]/F");
  tree->Branch("tk_qoverperr", &tk_qoverperr,"tk_qoverperr[tk_n]/F" );
  tree->Branch("tk_pterr", &tk_pterr,"tk_pterr[tk_n]/F" );
  tree->Branch("tk_etaerr", &tk_etaerr,"tk_etaerr[tk_n]/F" );
  tree->Branch("tk_phierr", &tk_phierr,"tk_phierr[tk_n]/F" );
  tree->Branch("tk_d0err", &tk_d0err,"tk_d0err[tk_n]/F" );
  tree->Branch("tk_dzerr", &tk_dzerr, "tk_dzerr[tk_n]/F");  
  tree->Branch("tk_hp_nvalid", &tk_hp_nvalid, "tk_hp_nvalid[tk_n]/I");  
  tree->Branch("tk_hp_nlost", &tk_hp_nlost, "tk_hp_nlost[tk_n]/I");  
  tree->Branch("tk_hp_nvalidpix", &tk_hp_nvalidpix, "tk_hp_nvalidpix[tk_n]/I");  
  tree->Branch("tk_hp_expin", &tk_hp_expin, "tk_hp_expin[tk_n]/I");  
  tree->Branch("tk_hp_expout", &tk_hp_expout, "tk_hp_expout[tk_n]/I"); 

  tree->Branch("tk_quality", &tk_quality, "tk_quality[tk_n]/I");
  tree->Branch("tk_algo", &tk_algo, "tk_algo[tk_n]/I");
}

//----------------------------------------------------------------------

bool GlobeTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  if (debug_level > 99)
    std::cout << "GlobeTracks: start " << std::endl;

  edm::Handle<reco::TrackCollection> tkH;
  //edm::Handle<reco::TrackCollection> tkH2;

  //if (nome == "ctf") // ONLY ONE TYPE
  iEvent.getByLabel(trackColl, tkH);
  //iEvent.getByLabel(trackColl2, tkH2);

  tk_p4->Clear();
  tk_vtx_pos->Clear();
  
  tk_n = 0;
  
  if (debug_level > 9)
    std::cout << "GlobeTracks: Track collection size: "<< tkH->size() << std::endl;
  //if (debug_level > 9)
  //  std::cout << "GlobeTracks: Track collection2 size: "<< tkH2->size() << std::endl;
 
  for(unsigned int i=0; i<tkH->size(); i++) {
    if (tk_n >= MAX_TRACKS) {
      std::cout << "GLobeTracks: WARNING TOO MANY TRACK: " << tkH->size() << " (allowed " << MAX_TRACKS << ")" << std::endl;
      break;
    }

    reco::TrackRef tk(tkH, i);
    // make the cuts
    if(gCUT->cut(*tk))
      continue; 
    // passed cuts
      

    new ((*tk_p4)[tk_n]) TLorentzVector();
    ((TLorentzVector *)tk_p4->At(tk_n))->SetXYZT(tk->px(), tk->py(), tk->pz(), tk->p());

    new ((*tk_vtx_pos)[tk_n]) TVector3();
    ((TVector3 *)tk_vtx_pos->At(tk_n))->SetXYZ(tk->vx(), tk->vy(), tk->vz());
    tk_charge[tk_n] = tk->charge();
    tk_nhits[tk_n] = tk->numberOfValidHits();
    tk_nlosthit[tk_n] = tk->numberOfLostHits() ;
    tk_dof[tk_n] = tk->ndof();
    tk_chi2[tk_n] = tk->chi2();
    tk_algo[tk_n] = (int)tk->algo();
    tk_quality[tk_n] = tk->qualityMask();
    
    tk_d0[tk_n] = tk->d0();
    tk_dz[tk_n] = tk->dz();
    tk_qoverperr[tk_n] = tk->qoverpError();
    tk_pterr[tk_n] = tk->ptError();
    tk_etaerr[tk_n] = tk->etaError();
    tk_phierr[tk_n] = tk->phiError();
    tk_dzerr[tk_n] = tk->dzError();
    tk_d0err[tk_n] = tk->d0Error();
    tk_cmsind[tk_n] = i;
    tk_tpind[tk_n] = -1;
    tk_hp_nvalid[tk_n] = tk->hitPattern().numberOfValidHits();
    tk_hp_nlost[tk_n] =  tk->hitPattern().numberOfLostHits();
    tk_hp_nvalidpix[tk_n] =  tk->hitPattern().numberOfValidPixelHits();

    //this is the same and it does not take forever
    tk_hp_expin[tk_n] = tk->trackerExpectedHitsInner().numberOfLostHits();
    tk_hp_expout[tk_n] = tk->trackerExpectedHitsOuter().numberOfLostHits();
    
    tk_n++;
  } // for i (loop over all tracks)

  //if (tk_n == 0)
  //  return false;
  
  return true;
}

//----------------------------------------------------------------------

//perform the match by associator (default is hits)
void GlobeTracks::GetAssociatedTrackingParticleIndex(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeTrackingParticles* tp) {

  edm::Handle<edm::View<reco::Track> > tkH;
  iEvent.getByLabel(trackColl, tkH);
  
  // get the TrackingParticles
  edm::Handle<TrackingParticleCollection> tpH;
  iEvent.getByLabel(tpColl, tpH);
  
  // get associator (by hits) for the track
  edm::ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get(assocLabel, theAssociator);
  
  // call the associator functions:
  reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(tkH, tpH, &iEvent);
  
  TrackingParticleRef associatedTrackingParticle;
  
  if(tp) for(int j=0;j<tk_n;j++) {
    edm::RefToBase<reco::Track> tk(tkH, (int)tk_cmsind[j]);

    if(recSimColl.find(tk) != recSimColl.end()) {

      // get the associated recoTrack      
      vector<pair<TrackingParticleRef, double> > simTracks = recSimColl[tk];
      associatedTrackingParticle = simTracks.begin()->first;
            
      // get the index of the associated recoTrack
      //for(TrackingParticleCollection::size_type i = 0; i < tpH->size(); ++i)
      for(int i = 0; i < tp->tp_n; ++i) {
        TrackingParticleRef trackingParticle(tpH, tp->tp_cmsind[i]);
        if( gCUT->cut(*trackingParticle) ) continue;
        
        if(associatedTrackingParticle == trackingParticle) {
          tk_tpind[j] = i;
          break;
        }
      }
    }
  }
}

//----------------------------------------------------------------------

std::pair<unsigned int, float> GlobeTracks::sharedHits(const reco::Track& trackA, const reco::Track& trackB) {
  
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

//----------------------------------------------------------------------
