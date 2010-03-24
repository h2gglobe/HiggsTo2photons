#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeVertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h" 

GlobeVertex::GlobeVertex(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  char a[100];
  sprintf(a, "VertexColl_%s", nome);
  vertexColl = iConfig.getParameter<edm::InputTag>(a);
  trackColl = iConfig.getParameter<edm::InputTag>("TrackColl");
  
  debug_level = iConfig.getParameter<int>("Debug_Level");
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
}

void GlobeVertex::defineBranch(TTree* tree) {

  vtx_xyz = new TClonesArray("TVector3", MAX_TRACKS);
  vtx_dxdydz = new TClonesArray("TVector3", MAX_TRACKS);
  vtx_vectorp3 = new TClonesArray("TVector3", MAX_TRACKS);
  
  char a1[50], a2[50];
  
  sprintf(a1, "vtx_%s_n", nome);
  sprintf(a2, "vtx_%s_n/I", nome);
  tree->Branch(a1, &vtx_n, a2);
  
  //CHECK add pttot, abspttot
  
  sprintf(a1, "vtx_%s_xyz", nome);
  tree->Branch(a1, "TClonesArray", &vtx_xyz, 32000, 0);
  
  sprintf(a1, "vtx_%s_dxdydz", nome);
  tree->Branch(a1, "TClonesArray", &vtx_dxdydz, 32000, 0);
  
  sprintf(a1, "vtx_%s_vectorp3", nome);
  tree->Branch(a1, "TClonesArray", &vtx_vectorp3, 32000, 0);
  
  sprintf(a1, "vtx_%s_x2dof", nome); 
  sprintf(a2, "vtx_%s_x2dof[vtx_%s_n]/F",nome,nome);
  tree->Branch(a1, &vtx_x2dof, a2);
  
  sprintf(a1, "vtx_%s_ndof", nome); 
  sprintf(a2, "vtx_%s_ndof[vtx_%s_n]/F",nome,nome);
  tree->Branch(a1, &vtx_ndof, a2);
  
  sprintf(a1, "vtx_%s_scalarpt", nome); 
  sprintf(a2, "vtx_%s_scalarpt[vtx_%s_n]/F",nome,nome);
  tree->Branch(a1, &vtx_scalarpt, a2);
  
  sprintf(a1, "vtx_%s_ntks", nome); 
  sprintf(a2, "vtx_%s_ntks[vtx_%s_n]/I",nome,nome);
  tree->Branch(a1, &vtx_ntks, a2);
  
  sprintf(a1, "vtx_%s_tkind", nome); 
  sprintf(a2, "vtx_%s_tkind[vtx_%s_n][%d]/I", nome, nome, MAX_VERTEX_TRACKS);
  tree->Branch(a1, &vtx_tkind, a2);
  
  sprintf(a1, "vtx_%s_tkweight", nome); 
  sprintf(a2, "vtx_%s_tkweight[vtx_%s_n][%d]/F", nome, nome, MAX_VERTEX_TRACKS);
  tree->Branch(a1, &vtx_tkweight, a2);
}

bool GlobeVertex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  if (debug_level > 99)
    std::cout << "GlobeVertex: start 1: "<< std::endl;
  
  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vertexColl, vtxH);
  
  if (debug_level > 99)
    std::cout << "GlobeVertex: start 2: "<< std::endl;
  
  edm::Handle<reco::TrackCollection> tkH;
  iEvent.getByLabel(trackColl, tkH);
  
  if (debug_level > 99)
    std::cout << "GlobeVertex: start 3: "<< std::endl;
  
  vtx_n = 0;
  
  vtx_xyz->Clear();
  vtx_dxdydz->Clear();
  vtx_vectorp3->Clear();
  
  if (debug_level > 9)
    std::cout << "GlobeVertex: Vertex collection size: "<< vtxH->size() << std::endl;
  
  for(unsigned int i=0; i<vtxH->size(); i++) {
    
    if (vtx_n >= MAX_VERTICES) {
      std::cout << "GlobeVertex: WARNING TOO MANY VERTEX: " << vtxH->size() << " (allowed " << MAX_VERTICES << ")" << std::endl;
      break;
    }
    
    reco::VertexRef vtx(vtxH, i);
    // Any cut on vertex?? gCUT ??  
    new((*vtx_xyz)[vtx_n]) TVector3();
    ((TVector3 *) vtx_xyz->At(vtx_n))->SetXYZ(vtx->x(), vtx->y(), vtx->z());
    
    new((*vtx_dxdydz)[vtx_n]) TVector3();
    ((TVector3 *)vtx_dxdydz->At(vtx_n))->SetXYZ(vtx->xError(), vtx->yError(), vtx->zError());
    
    math::XYZVector vecsumP3 = math::XYZVector(0.,0.,0.);
    vtx_scalarpt[vtx_n]=0.;
    
    for(std::vector<reco::TrackBaseRef>::const_iterator trkItr = vtx->tracks_begin();trkItr != vtx->tracks_end(); ++trkItr) {
      vtx_scalarpt[vtx_n] += (*trkItr)->pt();
      vecsumP3 += (*trkItr)->momentum();
    }
    new((*vtx_vectorp3)[vtx_n]) TVector3();
    ((TVector3 *)vtx_vectorp3->At(vtx_n))->SetXYZ(vecsumP3.x(), vecsumP3.y(), vecsumP3.z());
    
    
    if (strcmp(nome, "std") == 1) { 
      if (vtx->tracksSize() < MAX_VERTEX_TRACKS) {
        vtx_ntks[vtx_n] = vtx->tracksSize();
        
        if (debug_level > 9)
          std::cout << "GlobeVertex: Tracks size: "<< vtx->tracksSize() << std::endl;
        
        int k=0;
        std::vector<reco::TrackBaseRef>::const_iterator tk;
            
        for(tk=vtx->tracks_begin();tk!=vtx->tracks_end();++tk) {
          int index = -1;
          for(reco::TrackCollection::size_type j = 0; j<tkH->size(); ++j) {
            reco::TrackRef track(tkH, j);
            if(gCUT->cut(*track))continue; 
            if (&(**tk) == &(*track)) {
              //std::cout << "MWL: " << k << " " << index << " " << vtx->trackWeight(track) << std::endl;
              vtx_tkind[vtx_n][k] = index;
              vtx_tkweight[vtx_n][k] = vtx->trackWeight(track);
              break;
            }
            index++;
          }
          if(index == -1) {
            vtx_tkind[vtx_n][k] = -9999;
            vtx_tkweight[vtx_n][k] = -9999;
          }
          k++;
        }
        for(int k=vtx->tracksSize();k<MAX_VERTEX_TRACKS;++k) {
          vtx_tkweight[vtx_n][k] = -9999;
          vtx_tkind[vtx_n][k] = -9999;
        }
      } else { //end < MAX VERTEX TRACKS
        vtx_ntks[vtx_n] = 0;
        for(int k=0;k<MAX_VERTEX_TRACKS;++k) {
          vtx_tkweight[vtx_n][k] = -9999;
          vtx_tkind[vtx_n][k] = -9999;
        }
      }
    } else { //end "std"
      vtx_ntks[vtx_n] = 0;
      for(int i=0;i<MAX_VERTEX_TRACKS;++i) {
        vtx_tkweight[vtx_n][i] = -9999;
        vtx_tkind[vtx_n][i] = -9999;
      }
    }
    
    vtx_x2dof[vtx_n]=vtx->normalizedChi2();
    vtx_ndof[vtx_n]=vtx->ndof();
    vtx_n++;
  }
  
  if (debug_level > 99)
    std::cout << "GlobeVertex: start 1: "<< std::endl;
  
  return true;
}
