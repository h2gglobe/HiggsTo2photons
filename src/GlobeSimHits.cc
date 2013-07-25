#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimHits.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

GlobeSimHits::GlobeSimHits(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  simhitTIBLowColl			= iConfig.getParameter<edm::InputTag>("SimHitTIBLowColl");
  simhitTIDLowColl			= iConfig.getParameter<edm::InputTag>("SimHitTIDLowColl");
  simhitTECLowColl			= iConfig.getParameter<edm::InputTag>("SimHitTECLowColl");
  simhitTOBLowColl			= iConfig.getParameter<edm::InputTag>("SimHitTOBLowColl");
  simhitPixBarrelLowColl 	= iConfig.getParameter<edm::InputTag>("SimHitPixBarrelLowColl");
  simhitPixEndcapLowColl	= iConfig.getParameter<edm::InputTag>("SimHitPixEndcapLowColl");

  simtrackColl 	= iConfig.getParameter<edm::InputTag>("SimTrackColl");

  debug_level 		= iConfig.getParameter<int>("Debug_Level");

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig); 
}

void GlobeSimHits::defineBranch(GlobeAnalyzer* ana) {

  simhit_xyz= new TClonesArray("TVector3", MAX_SIMHITS);
  ana->Branch("simhit_n", &simhit_n, "simhit_n/I");
  ana->Branch("simhit_xyz", "TClonesArray", &simhit_xyz, 32000, 0);
  ana->Branch("simhit_pabs", &simhit_pabs, "simhit_pabs[simhit_n]/F");
  ana->Branch("simhit_eloss", &simhit_eloss, "simhit_eloss[simhit_n]/F");
  ana->Branch("simhit_subdet", &simhit_subdet, "simhit_subdet[simhit_n]/I");
  ana->Branch("simhit_pdgid", &simhit_pdgid, "simhit_pdgid[simhit_n]/I");
  ana->Branch("simhit_trkid", &simhit_trkid, "simhit_trkid[simhit_n]/I");
  ana->Branch("simhit_simtrkind", &simhit_simtrkind, "simhit_simtrkind[simhit_n]/I");
}

bool GlobeSimHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  // Tracker Si Sim Hits	 
  edm::Handle<edm::PSimHitContainer> simhitTIBLowH;
  iEvent.getByLabel(simhitTIBLowColl, simhitTIBLowH);
  edm::Handle<edm::PSimHitContainer> simhitTIDLowH;
  iEvent.getByLabel(simhitTIDLowColl, simhitTIDLowH);
  edm::Handle<edm::PSimHitContainer> simhitTOBLowH;
  iEvent.getByLabel(simhitTOBLowColl, simhitTOBLowH);
  edm::Handle<edm::PSimHitContainer> simhitTECLowH;
  iEvent.getByLabel(simhitTECLowColl, simhitTECLowH);
  // Tracker Pixel Sim Hits
  edm::Handle<edm::PSimHitContainer> simhitPixBarrelLowH;
  iEvent.getByLabel(simhitPixBarrelLowColl, simhitPixBarrelLowH);
  edm::Handle<edm::PSimHitContainer> simhitPixEndcapLowH;
  iEvent.getByLabel(simhitPixEndcapLowColl, simhitPixEndcapLowH);
	 
  theSimHits.clear();
  theSimHits.insert(theSimHits.end(),simhitTIBLowH->begin(),simhitTIBLowH->end());
  theSimHits.insert(theSimHits.end(),simhitTIDLowH->begin(),simhitTIDLowH->end());
  theSimHits.insert(theSimHits.end(),simhitTOBLowH->begin(),simhitTOBLowH->end());
  theSimHits.insert(theSimHits.end(),simhitTECLowH->begin(),simhitTECLowH->end());
  theSimHits.insert(theSimHits.end(),simhitPixBarrelLowH->begin(),simhitPixBarrelLowH->end());
  theSimHits.insert(theSimHits.end(),simhitPixEndcapLowH->begin(),simhitPixEndcapLowH->end());

  edm::Handle<edm::SimTrackContainer> stkH;
  iEvent.getByLabel(simtrackColl, stkH);
  theSimTracks.clear();
  theSimTracks.insert(theSimTracks.end(),stkH->begin(),stkH->end());

  SimTrack matchTrack;
  simhit_n = 0;
  simhit_xyz->Clear();
  for(std::vector<PSimHit>::iterator simhitItr = theSimHits.begin(); simhitItr != theSimHits.end(); ++simhitItr) {
    if(gCUT->cut(*simhitItr))continue;
    if (simhit_n >= MAX_SIMHITS) {
      if(MAX_SIMHITS>1) { //CHECK temporary because too slow and set it to 1
	std::cout << "GlobeSimHits: WARNING TOO MANY SIMHITS: " << theSimHits.size() << " (allowed " << MAX_SIMHITS << ")" << std::endl;
      }
      break;
    }
	 Int_t simtrk_n = 0;
	 Int_t simtrkindex = -1;
    for (std::vector<SimTrack>::iterator simtrackItr = theSimTracks.begin(); simtrackItr != theSimTracks.end(); ++simtrackItr){
      if(gCUT->cut(*simtrackItr))continue;
		if(simtrk_n >= MAX_SIMTRACKS) break;
	   if(simtrackItr->trackId() == simhitItr->trackId()) {
	     matchTrack = *simtrackItr;
		  simtrkindex = simtrk_n;
		  break;
		}
		simtrk_n++;
    }
	 Float_t lx = simhitItr->localPosition().x();
	 Float_t ly = simhitItr->localPosition().y();
	 Float_t lz = simhitItr->localPosition().z();
	 new ((*simhit_xyz)[simhit_n]) TVector3();
	 ((TVector3 *)simhit_xyz->At(simhit_n))->SetXYZ(lx, ly, lz);
	 simhit_pabs[simhit_n] = simhitItr->pabs();
	 simhit_eloss[simhit_n] = simhitItr->energyLoss();
	 DetId theDetUnitId(simhitItr->detUnitId());
	 simhit_subdet[simhit_n] = theDetUnitId.subdetId();
	 simhit_pdgid[simhit_n] = simhitItr->particleType();
	 simhit_trkid[simhit_n] = simhitItr->trackId();
	 simhit_simtrkind[simhit_n] = simtrkindex;
	 simhit_n++;
  }

  return true;

}

