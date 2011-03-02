#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTrackingParticles.h"

GlobeTrackingParticles::GlobeTrackingParticles(const edm::ParameterSet& iConfig, const char* n): nome(n)
{
    tpColl        = iConfig.getParameter<edm::InputTag>("TPColl");
    tvColl        = iConfig.getParameter<edm::InputTag>("TVColl");
    trackColl     = iConfig.getParameter<edm::InputTag>("TrackColl");
    assocLabel    = iConfig.getParameter<std::string>("AssocLabel");
    debug_level   = iConfig.getParameter<int>("Debug_Level");
    simtrackColl  = iConfig.getParameter<edm::InputTag>("SimTrackColl");
    simvertexColl     = iConfig.getParameter<edm::InputTag>("SimVertexColl");
    generatorColl = iConfig.getParameter<edm::InputTag>("GeneratorColl");

    // get cut thresholds
    gCUT = new GlobeCuts(iConfig); 
}

void GlobeTrackingParticles::defineBranch(TTree* tree) {

    // think about changing branch names for duplicate collections
    tp_p4 = new TClonesArray("TLorentzVector", MAX_TRACKINGPARTICLES);
    tp_vtx= new TClonesArray("TVector3", MAX_TRACKINGPARTICLES);
    tv_xyz= new TClonesArray("TVector3", MAX_TRACKINGPARTICLES);

    tree->Branch("tp_n", &tp_n, "tp_n/I");
    tree->Branch("tp_p4", "TClonesArray", &tp_p4, 32000, 0);
    tree->Branch("tp_vtx", "TClonesArray", &tp_vtx, 32000, 0);
    tree->Branch("tv_xyz", "TClonesArray", &tv_xyz, 32000, 0);
    tree->Branch("tp_pdgid", &tp_pdgid, "tp_pdgid[tp_n]/I");
    tree->Branch("tp_motherid", &tp_motherid, "tp_motherid[tp_n]/I");
    tree->Branch("tp_charge", &tp_charge, "tp_charge[tp_n]/I");
    tree->Branch("tp_tkind", &tp_tkind, "tp_tkind[tp_n]/I");
    tree->Branch("tp_genind", &tp_genind, "tp_genind[tp_n]/I");
    tree->Branch("tp_d0", &tp_d0, "tp_d0[tp_n]/D");
    tree->Branch("tp_dz", &tp_dz, "tp_dz[tp_n]/D");
}


bool GlobeTrackingParticles::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeTracks *tracks)
{
    using namespace std;

    // get the TrackingParticles
    edm::Handle<TrackingParticleCollection> trackingParticleHandle;
    iEvent.getByLabel(tpColl, trackingParticleHandle);

    // get the TrackingParticles
    edm::Handle<TrackingVertexCollection> trackingVertexHandle;
    iEvent.getByLabel(tvColl, trackingVertexHandle);

    // get the magnetic field 
    edm::ESHandle<MagneticField> mfH;
    iSetup.get<IdealMagneticFieldRecord>().get(mfH);
    const MagneticField *bField = mfH.product();

    // get the gen particles
    edm::Handle<edm::HepMCProduct> genEventHandle;
    iEvent.getByLabel(generatorColl, genEventHandle); 
   
    // get the reco::tracks
    edm::Handle<edm::View<reco::Track> > recoTrackHandle;
    iEvent.getByLabel(trackColl, recoTrackHandle);

    // get associator (by hits) for the track
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    iSetup.get<TrackAssociatorRecord>().get(assocLabel, theAssociator);

    // call the associator functions:
    reco::SimToRecoCollection simRecColl = theAssociator->associateSimToReco(recoTrackHandle,trackingParticleHandle, &iEvent);

    // Get the sim Tracks and vertices
    edm::Handle<edm::SimTrackContainer> stkH;
    iEvent.getByLabel(simtrackColl, stkH);
    theSimTracks = *stkH.product();

    edm::Handle<edm::SimVertexContainer> svtxH;
    iEvent.getByLabel(simvertexColl, svtxH);
    theSimVertices = *svtxH.product();

    ////////////////////////////////// Get the 1st tracking vertex /////////////////////////////////////
    TrackingVertexRef trackingVertex( trackingVertexHandle, 0);
    tv_xyz->Clear();
    new ((*tv_xyz)[0]) TVector3( trackingVertex->position().x(), trackingVertex->position().y(), trackingVertex->position().z());


    ////////////////////////////////// loop over tracking particles ////////////////////////////////////
    // initialize counter
    tp_n = 0;
    tp_vtx->Clear();
    tp_p4->Clear();
    for(TrackingParticleCollection::size_type i = 0; i < trackingParticleHandle->size(); ++i)
    {
        //initialization
        TrackingParticleRef trackingParticle(trackingParticleHandle, i);
        
        // apply cut and check for bounds
        if( gCUT->cut(*trackingParticle) ) continue;
        if(tp_n >= MAX_TRACKINGPARTICLES) 
        {
           if(MAX_TRACKINGPARTICLES>1)  //CHECK temporary because too slow and set it to 1
           { 
               std::cout << "GlobeTrackingParticles: WARNING TOO MANY TRACKINGPARTICLES: " 
                         << tp_n << " ( allowed " << MAX_TRACKINGPARTICLES << ")" << std::endl;
           }
           break;  // break out of for loop
        }
        
        // set tp variables
        new ((*tp_vtx)[tp_n]) TVector3( trackingParticle->vx(),
                                        trackingParticle->vy(), 
                                        trackingParticle->vz() );
        new ((*tp_p4)[tp_n]) TLorentzVector( trackingParticle->px(), 
                                             trackingParticle->py(), 
                                             trackingParticle->pz(), 
                                             trackingParticle->energy() );
        tp_pdgid[tp_n] = trackingParticle->pdgId();
        tp_motherid[tp_n] = GetMotherTrackingParticle(trackingParticle);
        tp_charge[tp_n] = trackingParticle->charge();
        tp_tkind[tp_n] = GetAssociatedRecoTrackIndex(trackingParticle, recoTrackHandle, simRecColl,tracks);

        //for(reco::TrackCollection::size_type i = 0; i < tkH->size(); ++i)
        tp_genind[tp_n] = GetAssociatedGenParticleIndex(trackingParticle, genEventHandle);
        tp_cmsind[tp_n] = i;
        CalcTrackingParticleImpactParam(tp_d0[tp_n], tp_dz[tp_n], trackingParticle, bField);

        // increment counter
        tp_n++;
    }
    
    return true;   
}    


Int_t GlobeTrackingParticles::GetAssociatedGenParticleIndex(TrackingParticleRef tp, const edm::Handle<edm::HepMCProduct>& genEvtH)
{
  
  int matchedGenParticleIndex = -1;
  HepMC::GenParticle gp;
  
  // assumption: one gen/sim particle per TP (true for CMSSW_1_6_X at least)
  if( tp->genParticle().size())
    {
      gp = **( tp->genParticle_begin() );
    }
  else
    {
      return -1;
    }
  
  const HepMC::GenEvent* myGenEvent = genEvtH->GetEvent();
  int genind = 0;
  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) 
    {
      if (genind >= MAX_GENERATOR) 
        {
          std::cout << "GlobeGenerator: WARNING TOO MANY Generator PARTICLES (allowed " << MAX_GENERATOR << ")" << std::endl;
          break;
        }
      if( gp == **it)
        {
          matchedGenParticleIndex = genind;
          break;
        }
      genind++;
    }
  
  return matchedGenParticleIndex;
}

//#ifndef CMSSW_VERSION_209_AND_210
////perform the match by associator (default is hits) 
//Int_t GlobeTrackingParticles::GetAssociatedRecoTrackIndex(TrackingParticleRef tp, const edm::Handle<reco::TrackCollection> &tkH, const reco::SimToRecoCollection &simRecColl, GlobeTracks *tracks)
//{
//  
//    reco::TrackRef associatedRecoTrack;
//    int matchedTrackIndex = -1;     // -1 unmatched        
//    
//    if(tracks) if(simRecColl.find(tp) != simRecColl.end())
//    {
//        // get the associated recoTrack
//        associatedRecoTrack = simRecColl[tp].begin()->first;
//    
//        // get the index of the associated recoTrack
//        for(int i = 0; i < tracks->tk_n; ++i)
//        {
//            reco::TrackRef recoTrack(tkH, tracks->tk_cmsind[i]);
//            if( gCUT->cut(*recoTrack) ) continue;
//
//            if(associatedRecoTrack == recoTrack)
//            {
//                matchedTrackIndex = i;
//                break;
//            }
//        } 
//    }
//  
//    return matchedTrackIndex;
//}
//#else
//perform the match by associator (default is hits) 
Int_t GlobeTrackingParticles::GetAssociatedRecoTrackIndex(TrackingParticleRef tp, const edm::Handle<edm::View<reco::Track> > &tkH, const reco::SimToRecoCollection &simRecColl, GlobeTracks *tracks) {
   
  edm::RefToBase<reco::Track> associatedRecoTrack;
  int matchedTrackIndex = -1;     // -1 unmatched        
  
  if(tracks) {
    if(simRecColl.find(tp) != simRecColl.end()) {
      
      // get the associated recoTrack
      associatedRecoTrack = simRecColl[tp].begin()->first;
      
      // get the index of the associated recoTrack
      for(int i = 0; i < tracks->tk_n; ++i) {
         edm::RefToBase<reco::Track> recoTrack(tkH, tracks->tk_cmsind[i]);
        if( gCUT->cut(*recoTrack) ) 
          continue;
        
        if(associatedRecoTrack == recoTrack) {
            matchedTrackIndex = i;
            break;
        }
      } 
    }
  }
 
  return matchedTrackIndex;
}
//#endif

void GlobeTrackingParticles::CalcTrackingParticleImpactParam(double &d0, double &dz, TrackingParticleRef tp, const MagneticField *bf) {
  try {
    GlobalPoint trackingParticleVertex( tp->vx(), tp->vy(), tp->vz() );
    GlobalVector trackingParticleP3(tp->px(), tp->py(), tp->pz() );
    TrackCharge trackingParticleCharge(tp->charge());
    
    FreeTrajectoryState ftsAtProduction( trackingParticleVertex, trackingParticleP3, trackingParticleCharge, bf );
    TSCPBuilderNoMaterial tscpBuilder;
    TrajectoryStateClosestToPoint tsAtClosestApproach =  tscpBuilder(ftsAtProduction,GlobalPoint(0,0,0));   //as in TrackProducerAlgorithm
    GlobalPoint v = tsAtClosestApproach.theState().position();
    GlobalVector p = tsAtClosestApproach.theState().momentum();
    
    
    d0 = -1.0 * ( (- v.x() * p.y() + v.y() * p.x() ) / p .perp() );
    dz = v.z() - ( v.x() * p.x() + v.y() * p.y() )/p.perp() * p.z()/p.perp();
  } catch (...) {
    std::cout << "ATTENTION GlobeTrackingParticles gives: TrajectoryStateClosestToPoint is invalid and cannot return any parameters!  "<<std::endl;
  }

  return;
}


Int_t GlobeTrackingParticles::GetMotherTrackingParticle(const TrackingParticle& tp)
{
    //std::cout << "GetMotherTrackingParticle: " << std::endl;
    if(tp.genParticle().size() > 0) 
        return GetMotherGenParticle( **tp.genParticle().begin() );
    else if(tp.genParticle().size() == 0 && tp.g4Tracks().size() > 0) 
        return GetMotherSimTrack( *tp.g4Tracks().begin() );
    else return -9999;
}



Int_t GlobeTrackingParticles::GetMotherGenParticle(const GenParticle& gp)
{
    //std::cout << "GetMotherGenParticle: " << std::endl;
    // assumption: only one mother of the gen particle

    for(HepMC::GenVertex::particles_in_const_iterator iter = gp.production_vertex()->particles_in_const_begin();
        iter != gp.production_vertex()->particles_in_const_end(); ++iter)
    {
        if( (*iter)->pdg_id() == gp.pdg_id() )
        {
            return GetMotherGenParticle(**iter);
        }
        else
        {
            return (*iter)->pdg_id();
        }
    }
    
    return -9999;
}


// This function is not general.  It assumes that there are no genTracks associated to the trackingParticle.  Also,
// I had to put limits on the index to prevent seg faults.  I'm not sure why I have to assume this.  
Int_t GlobeTrackingParticles::GetMotherSimTrack(const SimTrack &st) 
{   
    static int depth = 1;
    //std::cout << "GetMotherSimTrack: depth == " << depth << std::endl;
    // assumption: only one mother of the simTrack 

    SimVertex prodVtx;
    SimTrack mother;
   
    //cout << "vertIndex() == " << st.vertIndex() << " theSimVertices.size() == " << theSimVertices.size() << " depth == " << depth << endl;
    if(st.vertIndex() >= 0 && st.vertIndex() < (int)theSimVertices.size()) 
    {
        prodVtx = theSimVertices[st.vertIndex()];
        //cout << "parentIndex() == " << prodVtx.parentIndex() << " theSimTracks.size() == " << theSimTracks.size() << " depth == " << depth << endl;
        if(prodVtx.parentIndex() >= 0 && prodVtx.parentIndex() < (int)theSimTracks.size()) 
        {
            mother = theSimTracks[prodVtx.parentIndex()];
            //cout << "mother.type() == " << mother.type() << " st.type() == " << st.type() << " depth == " << depth << endl;
            if(mother.type() == st.type() && mother.trackId() != st.trackId()) {
                depth++;
                return GetMotherSimTrack(mother);
            }
            depth--;
            return mother.type();
        }
        else 
        {
            depth--;
            return -9999;
        }
    } 
    
    depth--;
    return -9999;
}

