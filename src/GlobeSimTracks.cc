#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

GlobeSimTracks::GlobeSimTracks(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  simtrackColl 	= iConfig.getParameter<edm::InputTag>("SimTrackColl");
  simvertexColl 	= iConfig.getParameter<edm::InputTag>("SimVertexColl");
  debug_level 		= iConfig.getParameter<int>("Debug_Level");

  track_and_vertex = iConfig.getParameter<bool>("doSimTrackPlusSimVertex");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig); 
}

void GlobeSimTracks::defineBranch(GlobeAnalyzer* ana) {

  // think about changing branch names for duplicate collections
  if (track_and_vertex) {
    simtrk_vtx= new TClonesArray("TVector3", MAX_SIMTRACKS);
    simtrk_p4= new TClonesArray("TLorentzVector", MAX_SIMTRACKS);
    
    ana->Branch("simtrk_n", &simtrk_n, "simtrk_n/I");
    ana->Branch("simtrk_p4", "TClonesArray", &simtrk_p4, 32000, 0);
    ana->Branch("simtrk_vtx", "TClonesArray", &simtrk_vtx, 32000, 0);
    ana->Branch("simtrk_pdgid", &simtrk_pdgid, "simtrk_pdgid[simtrk_n]/I");
    ana->Branch("simtrk_trkid", &simtrk_trkid, "simtrk_trkid[simtrk_n]/I");
    ana->Branch("simtrk_mothertrkid", &simtrk_mothertrkid, "simtrk_mothertrkid[simtrk_n]/I");
  }
  
  simvtx= new TClonesArray("TVector3", MAX_SIMTRACKS);
  ana->Branch("simvtx","TClonesArray", &simvtx, 32000, 0);
}

bool GlobeSimTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  SimVertex primVtx;
  //int npv=0;
  //int iPho=0;
  //int iPizero=0;
  std::vector<SimTrack> trkFromConversion;
  //stk_n = 0;

  edm::Handle<edm::SimTrackContainer> stkH;
  iEvent.getByLabel(simtrackColl, stkH);
  theSimTracks.clear();
  theSimTracks.insert(theSimTracks.end(),stkH->begin(),stkH->end());
  
  edm::Handle<edm::SimVertexContainer> svtxH;
  iEvent.getByLabel(simvertexColl, svtxH);
  theSimVertices.clear();
  theSimVertices.insert(theSimVertices.end(),svtxH->begin(),svtxH->end());

  fillGeantMap();

  //int iPV=-1; 
  std::vector<SimTrack>::iterator iFirstSimTk = theSimTracks.begin();
  if (  !(*iFirstSimTk).noVertex() ) {
    //iPV =  (*iFirstSimTk).vertIndex();
    int vtxId =   (*iFirstSimTk).vertIndex();
    primVtx = theSimVertices[vtxId];  
  }
  simvtx->Clear();

  if (  ( !( (*iFirstSimTk).noVertex() ) )) {
    if(svtxH->size()) {
      new ((*simvtx)[0]) TVector3(primVtx.position().x(), primVtx.position().y(), primVtx.position().z());
    }
    else {
      new ((*simvtx)[0]) TVector3(0., 0., 0.);
    }
  }
  else {
    new ((*simvtx)[0]) TVector3(0., 0., 0.);
  }

  if (track_and_vertex) {
    simtrk_n = 0;
    simtrk_vtx->Clear();
    simtrk_p4->Clear();
    for (std::vector<SimTrack>::iterator simtrackItr = theSimTracks.begin(); simtrackItr != theSimTracks.end(); ++simtrackItr){
      if(gCUT->cut(*simtrackItr))continue;
      if (simtrk_n >= MAX_SIMTRACKS) {
        if(MAX_SIMTRACKS>1) { //CHECK temporary because too slow and set it to 1
          std::cout << "GlobeSimTracks: WARNING TOO MANY SIMTRACKS: " << theSimTracks.size() << " (allowed " << MAX_SIMTRACKS << ")" << std::endl;
        }
        break;
      }
      
      float x = theSimVertices[simtrackItr->vertIndex()].position().x();
      float y = theSimVertices[simtrackItr->vertIndex()].position().y();
      float z = theSimVertices[simtrackItr->vertIndex()].position().z();
      new ((*simtrk_vtx)[simtrk_n]) TVector3(x, y, z);
      float px = simtrackItr->momentum().px();
      float py = simtrackItr->momentum().py();
      float pz = simtrackItr->momentum().pz();
      float en = simtrackItr->momentum().e();
      new ((*simtrk_p4)[simtrk_n]) TLorentzVector(px, py, pz, en);
      
      simtrk_pdgid[simtrk_n] = simtrackItr->type();
      simtrk_trkid[simtrk_n] = simtrackItr->trackId();
      SimTrack motherTrack = getMotherSimTrack(*simtrackItr);
      simtrk_mothertrkid[simtrk_n] = motherTrack.trackId();
      
      simtrk_n++;
    }
  }

  return true;
}



SimTrack GlobeSimTracks::getNextSimTrack(SimTrack &inputSimtrack) {

  Int_t simtrk_count=0;
  int inputVertexIndex= inputSimtrack.vertIndex();
  unsigned int inputTrackId= inputSimtrack.trackId();
  for (std::vector<SimTrack>::iterator candSimtrack = theSimTracks.begin(); candSimtrack != theSimTracks.end(); ++candSimtrack){
    if (gCUT->cut(*candSimtrack))continue;
    if (simtrk_count >= MAX_SIMTRACKS) { break; }
	 simtrk_count++;
    if ( candSimtrack->vertIndex() == inputVertexIndex )continue; 
    SimVertex vertex = theSimVertices[candSimtrack->vertIndex()];
    if ( vertex.parentIndex() != -1 ) {
      
      unsigned  motherGeantIndex = vertex.parentIndex(); 
      std::map<unsigned, unsigned >::iterator association = geantToIndex_.find( motherGeantIndex );
      int motherIndex=-1;
      if(association != geantToIndex_.end() )
	motherIndex = association->second;
      if ( theSimTracks[motherIndex].trackId() == inputTrackId )
	return *candSimtrack;
    }

  }

  return inputSimtrack;
}


SimTrack GlobeSimTracks::getSiblingSimTrack(SimTrack &inputSimtrack) {

  if(inputSimtrack.noVertex()) return inputSimtrack;
  int inputVertexIndex= inputSimtrack.vertIndex();
  SimVertex inputVertex = theSimVertices[inputVertexIndex];
  if(inputVertex.parentIndex() == -1)return inputSimtrack;
  unsigned int inputMotherGeantIndex = inputVertex.parentIndex(); 
  std::map<unsigned, unsigned >::iterator inputAssociation = geantToIndex_.find( inputMotherGeantIndex );
  if(inputAssociation == geantToIndex_.end() )return inputSimtrack;
  int inputMotherIndex = -1; 
  inputMotherIndex = inputAssociation->second;
  unsigned int inputMotherId= theSimTracks[inputMotherIndex].trackId();
  unsigned int inputTrackId= inputSimtrack.trackId();

  Int_t simtrk_count=0;
  for (std::vector<SimTrack>::iterator candSimtrack = theSimTracks.begin(); candSimtrack != theSimTracks.end(); ++candSimtrack){
    if ( gCUT->cut(*candSimtrack) )continue;
    if (simtrk_count >= MAX_SIMTRACKS) { break; }
	 simtrk_count++;
    if ( candSimtrack->vertIndex() != inputVertexIndex )continue; 
    if ( candSimtrack->trackId() == inputTrackId ) continue;

    SimVertex vertex = theSimVertices[candSimtrack->vertIndex()];
    if ( vertex.parentIndex() != -1 ) {

      unsigned  motherGeantIndex = vertex.parentIndex(); 
      std::map<unsigned, unsigned >::iterator association = geantToIndex_.find( motherGeantIndex );
      int motherIndex=-1;
      if(association != geantToIndex_.end() )
	     motherIndex = association->second;
      if ( theSimTracks[motherIndex].trackId() == inputMotherId )
	     return *candSimtrack;
    }

  }

  return inputSimtrack;
}

SimTrack GlobeSimTracks::getMotherSimTrack(SimTrack &inputSimtrack) {

  Int_t simtrk_count=0;
  unsigned int inputTrackId = inputSimtrack.trackId();
  for (std::vector<SimTrack>::iterator motherItr = theSimTracks.begin(); motherItr != theSimTracks.end(); ++motherItr){
    if ( gCUT->cut(*motherItr) )continue;
    if (simtrk_count >= MAX_SIMTRACKS) { break; }
	 simtrk_count++;
    SimTrack sameTrack = getNextSimTrack(*motherItr);
    if(sameTrack.trackId() == inputTrackId) return *motherItr;
  }
  return inputSimtrack;
}

void GlobeSimTracks::fillGeantMap() {
  //std::cout << "  GlobeSimTracks::fill " << std::endl;
  // Watch out there ! A SimVertex is in mm (stupid), 
  unsigned nVtx = theSimVertices.size();
  unsigned nTks = theSimTracks.size();
  // Empty event
  if ( nVtx == 0 ) return;
  // create a map associating geant particle id and position in the
  // event SimTrack vecto
  for( unsigned it=0; it<nTks; ++it ) {
    geantToIndex_[ theSimTracks[it].trackId() ] = it;
    //std::cout << " GlobeSimTracks::fill it " << it << " theSimTracks[it].trackId() " <<  theSimTracks[it].trackId() << std::endl;
  } 
}
