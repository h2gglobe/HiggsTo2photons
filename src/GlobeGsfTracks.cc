#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGsfTracks.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"

GlobeGsfTracks::GlobeGsfTracks(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  doAodSim     = iConfig.getParameter<bool>("doAodSim");
  storeGsfTracksOnlyIfElectrons     = iConfig.getParameter<bool>("storeGsfTracksOnlyIfElectrons");
  trackColl    = iConfig.getParameter<edm::InputTag>("GsfTrackColl");
  trackColl2   = iConfig.getParameter<edm::InputTag>("TrackColl");
  tpColl       = iConfig.getParameter<edm::InputTag>("TPColl");
  assocLabel   = iConfig.getParameter<std::string>("AssocLabel");
  debug_level  = iConfig.getParameter<int>("Debug_Level");
  electronColl = iConfig.getParameter<edm::InputTag>("ElectronColl_std");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeGsfTracks::defineBranch(TTree* tree) {

  
  gsf_tk_p4 = new TClonesArray("TLorentzVector", MAX_TRACKS);
  gsf_tk_vtx_pos = new TClonesArray("TVector3", MAX_TRACKS);
  gsf_tk_pinmode = new TClonesArray("TLorentzVector", MAX_TRACKS);
  gsf_tk_poutmode = new TClonesArray("TLorentzVector", MAX_TRACKS);

  tree->Branch("gsf_tk_n", &gsf_tk_n, "gsf_tk_n/I");
  tree->Branch("gsf_tk_p4", "TClonesArray", &gsf_tk_p4, 32000, 0);
  tree->Branch("gsf_tk_vtx_pos", "TClonesArray", &gsf_tk_vtx_pos, 32000, 0);
  tree->Branch("gsf_tk_nhits", &gsf_tk_nhits, "gsf_tk_nhits[gsf_tk_n]/I");
  tree->Branch("gsf_tk_charge", &gsf_tk_charge, "gsf_tk_charge[gsf_tk_n]/I");
  tree->Branch("gsf_tk_nlosthit", &gsf_tk_nlosthit,"gsf_tk_nlosthit[gsf_tk_n]/I" );
  tree->Branch("gsf_tk_tpind", &gsf_tk_tpind,"gsf_tk_tpind[gsf_tk_n]/I" );
  tree->Branch("gsf_tk_chi2", &gsf_tk_chi2,"gsf_tk_chi2[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_dof", &gsf_tk_dof, "gsf_tk_dof[gsf_tk_n]/F");
  tree->Branch("gsf_tk_d0", &gsf_tk_d0, "gsf_tk_d0[gsf_tk_n]/F");
  tree->Branch("gsf_tk_dz", &gsf_tk_dz, "gsf_tk_dz[gsf_tk_n]/F");
  tree->Branch("gsf_tk_qoverpinerr", &gsf_tk_qoverpinerr,"gsf_tk_qoverpinerr[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_qoverpouterr", &gsf_tk_qoverpouterr,"gsf_tk_qoverpouterr[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_pterr", &gsf_tk_pterr,"gsf_tk_pterr[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_etaerr", &gsf_tk_etaerr,"gsf_tk_etaerr[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_phierr", &gsf_tk_phierr,"gsf_tk_phierr[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_d0err", &gsf_tk_d0err,"gsf_tk_d0err[gsf_tk_n]/F" );
  tree->Branch("gsf_tk_dzerr", &gsf_tk_dzerr, "gsf_tk_dzerr[gsf_tk_n]/F");  
  tree->Branch("gsf_tk_hp_nvalid", &gsf_tk_hp_nvalid, "gsf_tk_hp_nvalid[gsf_tk_n]/I");  
  tree->Branch("gsf_tk_hp_nlost", &gsf_tk_hp_nlost, "gsf_tk_hp_nlost[gsf_tk_n]/I");  
  tree->Branch("gsf_tk_hp_nvalidpix", &gsf_tk_hp_nvalidpix, "gsf_tk_hp_nvalidpix[gsf_tk_n]/I");  
  tree->Branch("gsf_tk_hp_expin", &gsf_tk_hp_expin, "gsf_tk_hp_expin[gsf_tk_n]/I");  
  tree->Branch("gsf_tk_hp_expout", &gsf_tk_hp_expout, "gsf_tk_hp_expout[gsf_tk_n]/I");  
  tree->Branch("gsf_tk_pin", &gsf_tk_pin, "gsf_tk_pin[gsf_tk_n]/F");    
  tree->Branch("gsf_tk_pout", &gsf_tk_pout, "gsf_tk_pout[gsf_tk_n]/F");    
  tree->Branch("gsf_tk_fbrem", &gsf_tk_fbrem, "gsf_tk_fbrem[gsf_tk_n]/F");  
  tree->Branch("gsf_tk_pinmode", "TClonesArray", &gsf_tk_pinmode, 32000, 0);
  tree->Branch("gsf_tk_poutmode", "TClonesArray", &gsf_tk_poutmode, 32000, 0);

  tree->Branch("gsf_tk_tkind", &gsf_tk_tkind, "gsf_tk_tkind[gsf_tk_n]/I");
  tree->Branch("gsf_tk_shared", &gsf_tk_shared, "gsf_tk_shared[gsf_tk_n]/F");
}

bool GlobeGsfTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  
  if (debug_level > 99)
    std::cout << "GlobeGsfTracks: start " << std::endl;


  edm::Handle<reco::GsfElectronCollection> el;
  iEvent.getByLabel(electronColl, el);

  gsf_tk_p4->Clear();
  gsf_tk_vtx_pos->Clear();
  gsf_tk_pinmode->Clear();
  gsf_tk_poutmode->Clear();

  gsf_tk_n = 0;
  
  if (storeGsfTracksOnlyIfElectrons && el->size() < 1)
    return true;

  edm::Handle<reco::GsfTrackCollection> tkH;
  iEvent.getByLabel(trackColl, tkH);
  edm::Handle<reco::TrackCollection> tkH2;
  iEvent.getByLabel(trackColl2, tkH2);

  if (debug_level > 9)
    std::cout << "GlobeGsfTracks: Track collection size: "<< tkH->size() << std::endl;
 
  for(unsigned int i=0; i<tkH->size(); i++) {
    if (gsf_tk_n >= MAX_TRACKS) {
      std::cout << "GLobeGsfTracks: WARNING TOO MANY TRACK: " << tkH->size() << " (allowed " << MAX_TRACKS << ")" << std::endl;
      break;
    }

    reco::GsfTrackRef tk(tkH, i);
      
    new ((*gsf_tk_p4)[gsf_tk_n]) TLorentzVector();
    ((TLorentzVector *)gsf_tk_p4->At(gsf_tk_n))->SetXYZT(tk->px(), tk->py(), tk->pz(), tk->p());

    new ((*gsf_tk_vtx_pos)[gsf_tk_n]) TVector3();
    ((TVector3 *)gsf_tk_vtx_pos->At(gsf_tk_n))->SetXYZ(tk->vx(), tk->vy(), tk->vz());
    gsf_tk_charge[gsf_tk_n] = tk->charge();
    gsf_tk_nhits[gsf_tk_n] = tk->numberOfValidHits();
    gsf_tk_nlosthit[gsf_tk_n] = tk->numberOfLostHits() ;
    gsf_tk_dof[gsf_tk_n] = tk->ndof();
    gsf_tk_chi2[gsf_tk_n] = tk->chi2();
    
    gsf_tk_d0[gsf_tk_n] = tk->d0();
    gsf_tk_dz[gsf_tk_n] = tk->dz();

    gsf_tk_pterr[gsf_tk_n] = tk->ptError();
    gsf_tk_etaerr[gsf_tk_n] = tk->etaError();
    gsf_tk_phierr[gsf_tk_n] = tk->phiError();
    gsf_tk_dzerr[gsf_tk_n] = tk->dzError();
    gsf_tk_d0err[gsf_tk_n] = tk->d0Error();
    gsf_tk_cmsind[gsf_tk_n] = i;
    gsf_tk_tpind[gsf_tk_n] = -1;
    gsf_tk_hp_nvalid[gsf_tk_n] = tk->hitPattern().numberOfValidHits();
    gsf_tk_hp_nlost[gsf_tk_n] =  tk->hitPattern().numberOfLostHits();
    gsf_tk_hp_nvalidpix[gsf_tk_n] =  tk->hitPattern().numberOfValidPixelHits();

    if (!doAodSim) {
      if (cacheIDMagField_!=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier()){
	
	cacheIDMagField_=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
	iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
      }
      
      if (cacheIDTDGeom_!=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier()){
	cacheIDTDGeom_=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier();
	iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
      }

      const MultiTrajectoryStateTransform *mtsTransform_ = new MultiTrajectoryStateTransform(trackerHandle_.product(),theMagField.product());
      const MultiTrajectoryStateMode *mtsMode_ = new MultiTrajectoryStateMode();
      TrajectoryStateOnSurface innTSOS_;
      TrajectoryStateOnSurface outTSOS_;
      
      innTSOS_ = mtsTransform_->innerStateOnSurface(*tk);
      outTSOS_ = mtsTransform_->outerStateOnSurface(*tk);
      GlobalVector innMom, outMom;
      mtsMode_->momentumFromModeCartesian(innTSOS_, innMom);
      mtsMode_->momentumFromModeCartesian(outTSOS_,outMom);
      gsf_tk_pin[gsf_tk_n] = innMom.mag();
      gsf_tk_pout[gsf_tk_n] = outMom.mag();
      gsf_tk_fbrem[gsf_tk_n] = (innMom.mag() - outMom.mag())/innMom.mag();
      
      MultiGaussianState1D qpstate = MultiGaussianStateTransform::innerMultiState1D(*tk, 0);
      GaussianSumUtilities1D qputils(qpstate);
      double qpvar_in = qputils.mode().variance();
      
      MultiGaussianState1D qpstate_out = MultiGaussianStateTransform::outerMultiState1D(*tk, 0);
      GaussianSumUtilities1D qputils_out(qpstate_out);
      double qpvar_out = qputils_out.mode().variance();
      gsf_tk_qoverpinerr[gsf_tk_n] =  qpvar_in; 
      gsf_tk_qoverpouterr[gsf_tk_n] = qpvar_out;

      delete mtsTransform_;
      delete mtsMode_;
    }

    gsf_tk_hp_expin[gsf_tk_n] = tk->trackerExpectedHitsInner().numberOfHits();
    gsf_tk_hp_expout[gsf_tk_n] = tk->trackerExpectedHitsOuter().numberOfHits();
    
    if (!doAodSim) {
      std::pair<reco::TrackRef,float> ctfRef = getCtfTrackRef(tk, tkH2);
      if (ctfRef.second > 0.) {
	gsf_tk_tkind[gsf_tk_n] = ctfRef.first.index();
	gsf_tk_shared[gsf_tk_n] = ctfRef.second;
      } else {
	gsf_tk_tkind[gsf_tk_n] = -1;
	gsf_tk_shared[gsf_tk_n] = -1;
      }
    }      
    
    gsf_tk_n++;
  }

  return true;
}

//perform the match by associator (default is hits)
void GlobeGsfTracks::GetAssociatedTrackingParticleIndex(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeTrackingParticles* tp) {

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
  
  if(tp) for(int j=0;j<gsf_tk_n;j++) {
    edm::RefToBase<reco::Track> tk(tkH, (int)gsf_tk_cmsind[j]);

    if(recSimColl.find(tk) != recSimColl.end()) {

      // get the associated recoTrack      
      vector<pair<TrackingParticleRef, double> > simTracks = recSimColl[tk];
      associatedTrackingParticle = simTracks.begin()->first;
            
      // get the index of the associated recoTrack
      for(int i = 0; i < tp->tp_n; ++i) {
        TrackingParticleRef trackingParticle(tpH, tp->tp_cmsind[i]);
        if( gCUT->cut(*trackingParticle) ) continue;
        
        if(associatedTrackingParticle == trackingParticle) {
          gsf_tk_tpind[j] = i;
          break;
        }
      }
    }
  }
}

std::pair<reco::TrackRef,float> GlobeGsfTracks::getCtfTrackRef
 ( const reco::GsfTrackRef& gsfTrackRef, edm::Handle<reco::TrackCollection> ctfTracksH )
 {
  float maxFracShared = 0;
  reco::TrackRef ctfTrackRef = reco::TrackRef() ;
  const reco::TrackCollection * ctfTrackCollection = ctfTracksH.product() ;

  // get the Hit Pattern for the gsfTrack
  const reco::HitPattern& gsfHitPattern = gsfTrackRef->hitPattern();

  unsigned int counter ;
  reco::TrackCollection::const_iterator ctfTkIter ;
  for ( ctfTkIter = ctfTrackCollection->begin() , counter = 0 ;
        ctfTkIter != ctfTrackCollection->end() ; ctfTkIter++, counter++ )
   {

    double dEta = gsfTrackRef->eta() - ctfTkIter->eta();
    double dPhi = gsfTrackRef->phi() - ctfTkIter->phi();
    double pi = acos(-1.);
    if(fabs(dPhi) > pi) dPhi = 2*pi - fabs(dPhi);

    // dont want to look at every single track in the event!
    if(sqrt(dEta*dEta + dPhi*dPhi) > 0.3) continue;

    unsigned int shared = 0 ;
    int gsfHitCounter = 0 ;
    int numGsfInnerHits = 0 ;
    int numCtfInnerHits = 0 ;
    // get the CTF Track Hit Pattern
    const reco::HitPattern& ctfHitPattern = ctfTkIter->hitPattern() ;

    trackingRecHit_iterator elHitsIt ;
    for ( elHitsIt = gsfTrackRef->recHitsBegin() ;
          elHitsIt != gsfTrackRef->recHitsEnd() ;
          elHitsIt++, gsfHitCounter++ )
     {
      if(!((**elHitsIt).isValid()))  //count only valid Hits
       { continue ; }

      // look only in the pixels/TIB/TID
      uint32_t gsfHit = gsfHitPattern.getHitPattern(gsfHitCounter) ;
      if (!(gsfHitPattern.pixelHitFilter(gsfHit) ||
          gsfHitPattern.stripTIBHitFilter(gsfHit) ||
          gsfHitPattern.stripTIDHitFilter(gsfHit) ) )
       { continue ; }

      numGsfInnerHits++ ;

      int ctfHitsCounter = 0 ;
      numCtfInnerHits = 0 ;
      trackingRecHit_iterator ctfHitsIt ;
      for ( ctfHitsIt = ctfTkIter->recHitsBegin() ;
            ctfHitsIt != ctfTkIter->recHitsEnd() ;
            ctfHitsIt++, ctfHitsCounter++ )
       {
        if(!((**ctfHitsIt).isValid())) //count only valid Hits!
         { continue ; }

      uint32_t ctfHit = ctfHitPattern.getHitPattern(ctfHitsCounter);
      if( !(ctfHitPattern.pixelHitFilter(ctfHit) ||
            ctfHitPattern.stripTIBHitFilter(ctfHit) ||
            ctfHitPattern.stripTIDHitFilter(ctfHit) ) )
       { continue ; }

      numCtfInnerHits++ ;

        if( (**elHitsIt).sharesInput(&(**ctfHitsIt),TrackingRecHit::all) )
         {
          shared++ ;
          break ;
         }

       } //ctfHits iterator

     } //gsfHits iterator

    if ((numGsfInnerHits==0)||(numCtfInnerHits==0))
     { continue ; }

    if ( static_cast<float>(shared)/std::min(numGsfInnerHits,numCtfInnerHits) > maxFracShared )
     {
      maxFracShared = static_cast<float>(shared)/std::min(numGsfInnerHits, numCtfInnerHits);
      ctfTrackRef = reco::TrackRef(ctfTracksH,counter);
     }

   } //ctfTrack iterator

  return make_pair(ctfTrackRef,maxFracShared) ;
 }
