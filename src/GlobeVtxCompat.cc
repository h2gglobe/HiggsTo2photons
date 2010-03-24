#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeVtxCompat.h"

GlobeVtxCompat::GlobeVtxCompat(const edm::ParameterSet& iConfig) {

    debug_level = iConfig.getParameter<int>("Debug_Level");
    trackColl =  iConfig.getParameter<edm::InputTag>("TrackColl");

    // get cut thresholds
    gCUT = new GlobeCuts(iConfig);
}

void GlobeVtxCompat::defineBranch(TTree* tree) {

  tree->Branch("vtxcomp_n", &vtxcomp_n,"vtxcomp_n/I");
  tree->Branch("vtxcomp_lepn", &vtxcomp_lepn,"vtxcomp_lepn[vtxcomp_n]/I");
  tree->Branch("vtxcomp_pdgid", &vtxcomp_pdgid,"vtxcomp_pdgid[vtxcomp_n][4]/I");
  tree->Branch("vtxcomp_lepind", &vtxcomp_lepind,"vtxcomp_lepind[vtxcomp_n][4]/I");
  tree->Branch("vtxcomp_klmn_isvalid",  &vtxcomp_klmn_isvalid  ,"vtxcomp_klmn_isvalid[vtxcomp_n]/O");
  tree->Branch("vtxcomp_klmn_x2dof",    &vtxcomp_klmn_x2dof    ,"vtxcomp_klmn_x2dof[vtxcomp_n]/F");
  tree->Branch("vtxcomp_klmn_ndof",     &vtxcomp_klmn_ndof     ,"vtxcomp_klmn_ndof[vtxcomp_n]/I");
  tree->Branch("vtxcomp_klmn_x2prob",   &vtxcomp_klmn_x2prob   ,"vtxcomp_klmn_x2prob[vtxcomp_n]/F");
  tree->Branch("vtxcomp_klmn_tk_chi2",  &vtxcomp_klmn_tk_chi2  ,"vtxcomp_klmn_tk_chi2[vtxcomp_n][4]/F");
  //tree->Branch("vtxcomp_adpt_isvalid",  &vtxcomp_adpt_isvalid  ,"vtxcomp_adpt_isvalid[vtxcomp_n]/O");
  //tree->Branch("vtxcomp_adpt_x2dof",    &vtxcomp_adpt_x2dof    ,"vtxcomp_adpt_x2dof[vtxcomp_n]/F");
  //tree->Branch("vtxcomp_adpt_ndof",     &vtxcomp_adpt_ndof     ,"vtxcomp_adpt_ndof[vtxcomp_n]/I");
  //tree->Branch("vtxcomp_adpt_x2prob",   &vtxcomp_adpt_x2prob   ,"vtxcomp_adpt_x2prob[vtxcomp_n]/F");
  //tree->Branch("vtxcomp_adpt_tk_chi2",  &vtxcomp_adpt_tk_chi2  ,"vtxcomp_adpt_tk_chi2[vtxcomp_n][4]/F");
  //tree->Branch("vtxcomp_adpt_tk_weight",&vtxcomp_adpt_tk_weight,"vtxcomp_adpt_tk_weight[vtxcomp_n][4]/F");
}

bool GlobeVtxCompat::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeLeptons *lep, GlobeElectrons *el, GlobeMuons *mu, GlobeTracks *tk) {

    if(!lep || !el || !mu || !tk) {
        std::cout << "Collection Missing to Do VtxCompat" << std::endl;
        return 0;
    }

    vtxcomp_n = 0;
    
    edm::ESHandle<MagneticField> B;
    iSetup.get<IdealMagneticFieldRecord>().get(B);

    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByLabel(trackColl, tkH);

    std::vector <reco::TransientTrack> vecOfTrks;

    KalmanVertexFitter fitterKalman(true); //track refit
    TransientVertex myVertexKalman;

//#ifdef CMSSW_VERSION_168
//    KalmanVertexTrackCompatibilityEstimator trkVtxFitKalman;
//#endif

//#ifdef CMSSW_VERSION_180_AND_209
    KalmanVertexTrackCompatibilityEstimator<5> trkVtxFitKalman;
//#endif
    
    // Already commented
    //AdaptiveVertexFitter fitterAdaptive; 
    //TransientVertex myVertexAdaptive;
    //KalmanVertexTrackCompatibilityEstimator trkVtxFitAdaptive;

    lepton lepTemp;
    std::vector<lepton> leps;

    for(int j=0;j<lep->lpt_emu_n;++j) { //Begin Lepton Loop List
      if(abs(lep->lpt_pdgid[j]) != 11 && abs(lep->lpt_pdgid[j]) != 13 ) {
        std::cout << "Error, Expected Lepton is not Electron or Muon: " << lep->lpt_pdgid[j] << std::endl;
        continue;
      }
      
      if( abs(lep->lpt_pdgid[j]) == 11 )
        lepTemp.trkind = tk->tk_cmsind[el->el_tkind[lep->lpt_ind[j]]];
      else if( abs(lep->lpt_pdgid[j]) == 13 )
        lepTemp.trkind = tk->tk_cmsind[mu->mu_tkind[lep->lpt_ind[j]]];
      lepTemp.lepind = j;
      leps.push_back(lepTemp);
    }
    
    for(unsigned int i=0;i<leps.size();++i) {
      for(unsigned int j=i+1;j<leps.size();++j) {

        if (tkH->size() < 1) {
          std::cout << "No Tracks in the event...:-((" << std::endl;
          return 0;
        }
        
        if(vtxcomp_n > MAX_VTXCOMP-1) {
          std::cout << "Holy Crap, Too many vtxcompat: " << std::endl;
          return 0;
        }
        
        try {
          //Do All combinations of 2:
          vecOfTrks.clear();
          reco::TrackRef t1(tkH, i); 
          vecOfTrks.push_back(reco::TransientTrack(t1,B.product()));
          
          reco::TrackRef t2(tkH, j); 
          vecOfTrks.push_back(reco::TransientTrack(t2,B.product()));
          
          //myVertexAdaptive = fitterAdaptive.vertex(vecOfTrks);
          myVertexKalman = fitterKalman.vertex(vecOfTrks);

          vtxcomp_lepn[vtxcomp_n]         = 2;
          vtxcomp_pdgid[vtxcomp_n][0]     = lep->lpt_pdgid[leps[i].lepind];
          vtxcomp_lepind[vtxcomp_n][0]    = leps[i].lepind;
          vtxcomp_pdgid[vtxcomp_n][1]     = lep->lpt_pdgid[leps[j].lepind];
          vtxcomp_lepind[vtxcomp_n][1]    = leps[j].lepind;
          vtxcomp_pdgid[vtxcomp_n][2]     = -9999;
          vtxcomp_lepind[vtxcomp_n][2]    = -1;
          vtxcomp_pdgid[vtxcomp_n][3]     = -9999;
          vtxcomp_lepind[vtxcomp_n][3]    = -1;
        } catch(...) {
          std::cout << "Tk Pt=0" << std::endl;
        }
        if( myVertexKalman.isValid() ) {
          vtxcomp_klmn_isvalid[vtxcomp_n] = myVertexKalman.isValid(); 
          vtxcomp_klmn_x2dof[vtxcomp_n] = myVertexKalman.totalChiSquared() / myVertexKalman.degreesOfFreedom();
          vtxcomp_klmn_ndof[vtxcomp_n] = (int)myVertexKalman.degreesOfFreedom();
          vtxcomp_klmn_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_klmn_x2dof[vtxcomp_n],vtxcomp_klmn_ndof[vtxcomp_n]);
#ifdef CMSSW_VERSION_210
          vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]);
          vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]);
          vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
#else
          vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]).second;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]).second;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
#endif          

        } else {
          vtxcomp_klmn_isvalid[vtxcomp_n] = 0;
          vtxcomp_klmn_x2dof[vtxcomp_n] = -9999;
          vtxcomp_klmn_ndof[vtxcomp_n] = -9999;
          vtxcomp_klmn_x2prob[vtxcomp_n] = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = -9999;
          vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
        }
        /*
        if( myVertexAdaptive.isValid() ) {
        vtxcomp_adpt_isvalid[vtxcomp_n] = myVertexAdaptive.isValid(); 
        vtxcomp_adpt_x2dof[vtxcomp_n] = myVertexAdaptive.totalChiSquared() / myVertexAdaptive.degreesOfFreedom();
        vtxcomp_adpt_ndof[vtxcomp_n] = (int)myVertexAdaptive.degreesOfFreedom();
        vtxcomp_adpt_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_adpt_x2dof[vtxcomp_n],vtxcomp_adpt_ndof[vtxcomp_n]);
        vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[0]);
        vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[1]);
        vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = -9999;
          vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = myVertexAdaptive.trackWeight(vecOfTrks[0]);
          vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = myVertexAdaptive.trackWeight(vecOfTrks[1]);
          vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = -9999;
        } else {
          vtxcomp_adpt_isvalid[vtxcomp_n] = 0;
          vtxcomp_adpt_x2dof[vtxcomp_n] = -9999;
          vtxcomp_adpt_ndof[vtxcomp_n] = -9999;
          vtxcomp_adpt_x2prob[vtxcomp_n] = -9999;
          vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = -9999;
          vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = -9999;
          vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = -9999;
          vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = -9999;
          vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = -9999;
        }
        */
        vtxcomp_n++;
        for(unsigned int k=j+1;k<leps.size();++k) {
          
          if(vtxcomp_n > MAX_VTXCOMP-1) {
            std::cout << "Holy Crap, Too many vtxcompat: " << std::endl;
            return 0;
          }
          
          //Do All combinations of 3:
          vecOfTrks.clear();
          reco::TrackRef t1(tkH, i); 
          vecOfTrks.push_back(reco::TransientTrack(t1,B.product()));
          
          reco::TrackRef t2(tkH, j); 
          vecOfTrks.push_back(reco::TransientTrack(t2,B.product()));
          
          reco::TrackRef t3(tkH, k); 
          vecOfTrks.push_back(reco::TransientTrack(t3,B.product()));
          //myVertexAdaptive = fitterAdaptive.vertex(vecOfTrks);
          myVertexKalman = fitterKalman.vertex(vecOfTrks);
          
          vtxcomp_lepn[vtxcomp_n]         = 3;
          vtxcomp_pdgid[vtxcomp_n][0]     = lep->lpt_pdgid[leps[i].lepind];
          vtxcomp_lepind[vtxcomp_n][0]    = leps[i].lepind;
          vtxcomp_pdgid[vtxcomp_n][1]     = lep->lpt_pdgid[leps[j].lepind];
          vtxcomp_lepind[vtxcomp_n][1]    = leps[j].lepind;
          vtxcomp_pdgid[vtxcomp_n][2]     = lep->lpt_pdgid[leps[k].lepind];
          vtxcomp_lepind[vtxcomp_n][2]    = leps[k].lepind;
          vtxcomp_pdgid[vtxcomp_n][3]     = -9999;
          vtxcomp_lepind[vtxcomp_n][3]    = -1;
	  //std::cout << "VTX6" << std::endl;
          if( myVertexKalman.isValid() ) {
            vtxcomp_klmn_isvalid[vtxcomp_n] = myVertexKalman.isValid(); 
            vtxcomp_klmn_x2dof[vtxcomp_n] = myVertexKalman.totalChiSquared() / myVertexKalman.degreesOfFreedom();
            vtxcomp_klmn_ndof[vtxcomp_n] = (int)myVertexKalman.degreesOfFreedom();
            vtxcomp_klmn_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_klmn_x2dof[vtxcomp_n],vtxcomp_klmn_ndof[vtxcomp_n]);
#ifdef CMSSW_VERSION_210
            vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]);
            vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]);
            vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[2]);
            vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
#else
            vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]).second;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]).second;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[2]).second;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
#endif
          } else {
            vtxcomp_klmn_isvalid[vtxcomp_n] = 0;
            vtxcomp_klmn_x2dof[vtxcomp_n] = -9999;
            vtxcomp_klmn_ndof[vtxcomp_n] = -9999;
            vtxcomp_klmn_x2prob[vtxcomp_n] = -9999;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = -9999;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = -9999;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = -9999;
            vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
          }
          
                /*if( myVertexAdaptive.isValid() ) {
                    vtxcomp_adpt_isvalid[vtxcomp_n] = myVertexAdaptive.isValid(); 
                    vtxcomp_adpt_x2dof[vtxcomp_n] = myVertexAdaptive.totalChiSquared() / myVertexAdaptive.degreesOfFreedom();
                    vtxcomp_adpt_ndof[vtxcomp_n] = (int)myVertexAdaptive.degreesOfFreedom();
                    vtxcomp_adpt_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_adpt_x2dof[vtxcomp_n],vtxcomp_adpt_ndof[vtxcomp_n]);
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[0]);
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[1]);
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[2]);
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = -9999;
                    vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = myVertexAdaptive.trackWeight(vecOfTrks[0]);
                    vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = myVertexAdaptive.trackWeight(vecOfTrks[1]);
                    vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = myVertexAdaptive.trackWeight(vecOfTrks[2]);
                    vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = -9999;
                } else {
                    vtxcomp_adpt_isvalid[vtxcomp_n] = 0;
                    vtxcomp_adpt_x2dof[vtxcomp_n] = -9999;
                    vtxcomp_adpt_ndof[vtxcomp_n] = -9999;
                    vtxcomp_adpt_x2prob[vtxcomp_n] = -9999;
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = -9999;
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = -9999;
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = -9999;
                    vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = -9999;
                    vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = -9999;
                    vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = -9999;
                    vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = -9999;
                    vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = -9999;
                }*/
          vtxcomp_n++;
          for(unsigned int l=k+1;l<leps.size();++l) {
            
            if(vtxcomp_n > MAX_VTXCOMP-1) {
              std::cout << "Holy Crap, Too many vtxcompat: " << std::endl;
              return 0;
            }
            
            //Do All combinations of 4:
            vecOfTrks.clear();
            reco::TrackRef t1(tkH, i); 
            vecOfTrks.push_back(reco::TransientTrack(t1,B.product()));
            
            reco::TrackRef t2(tkH, j); 
            vecOfTrks.push_back(reco::TransientTrack(t2,B.product()));
            
            reco::TrackRef t3(tkH, k); 
            vecOfTrks.push_back(reco::TransientTrack(t3,B.product()));
            
            reco::TrackRef t4(tkH, l); 
            vecOfTrks.push_back(reco::TransientTrack(t4,B.product()));
            
            //myVertexAdaptive = fitterAdaptive.vertex(vecOfTrks);
            myVertexKalman = fitterKalman.vertex(vecOfTrks);

            vtxcomp_lepn[vtxcomp_n]         = 4;
            vtxcomp_pdgid[vtxcomp_n][0]     = lep->lpt_pdgid[leps[i].lepind];
            vtxcomp_lepind[vtxcomp_n][0]    = leps[i].lepind;
            vtxcomp_pdgid[vtxcomp_n][1]     = lep->lpt_pdgid[leps[j].lepind];
            vtxcomp_lepind[vtxcomp_n][1]    = leps[j].lepind;
            vtxcomp_pdgid[vtxcomp_n][2]     = lep->lpt_pdgid[leps[k].lepind];
            vtxcomp_lepind[vtxcomp_n][2]    = leps[k].lepind;
            vtxcomp_pdgid[vtxcomp_n][3]     = lep->lpt_pdgid[leps[l].lepind];
            vtxcomp_lepind[vtxcomp_n][3]    = leps[l].lepind;
            
            if( myVertexKalman.isValid() ) {
              vtxcomp_klmn_isvalid[vtxcomp_n] = myVertexKalman.isValid(); 
              vtxcomp_klmn_x2dof[vtxcomp_n] = myVertexKalman.totalChiSquared() / myVertexKalman.degreesOfFreedom();
              vtxcomp_klmn_ndof[vtxcomp_n] = (int)myVertexKalman.degreesOfFreedom();
              vtxcomp_klmn_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_klmn_x2dof[vtxcomp_n],vtxcomp_klmn_ndof[vtxcomp_n]);
#ifdef CMSSW_VERSION_210
              vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]);
              vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]);
              vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[2]);
              vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[3]);
#else
              vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[0]).second;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[1]).second;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[2]).second;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = trkVtxFitKalman.estimate(myVertexKalman,vecOfTrks[3]).second;
#endif
            } else {
              vtxcomp_klmn_isvalid[vtxcomp_n] = 0;
              vtxcomp_klmn_x2dof[vtxcomp_n] = -9999;
              vtxcomp_klmn_ndof[vtxcomp_n] = -9999;
              vtxcomp_klmn_x2prob[vtxcomp_n] = -9999;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][0]  = -9999;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][1]  = -9999;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][2]  = -9999;
              vtxcomp_klmn_tk_chi2[vtxcomp_n][3]  = -9999;
            }
        
                    /*if( myVertexAdaptive.isValid() ) {
                        vtxcomp_adpt_isvalid[vtxcomp_n] = myVertexAdaptive.isValid(); 
                        vtxcomp_adpt_x2dof[vtxcomp_n] = myVertexAdaptive.totalChiSquared() / myVertexAdaptive.degreesOfFreedom();
                        vtxcomp_adpt_ndof[vtxcomp_n] = (int)myVertexAdaptive.degreesOfFreedom();
                        vtxcomp_adpt_x2prob[vtxcomp_n] = TMath::Prob(vtxcomp_adpt_x2dof[vtxcomp_n],vtxcomp_adpt_ndof[vtxcomp_n]);
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[0]);
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[1]);
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[2]);
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = trkVtxFitAdaptive.estimate(myVertexAdaptive,vecOfTrks[3]);
                        vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = myVertexAdaptive.trackWeight(vecOfTrks[0]);
                        vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = myVertexAdaptive.trackWeight(vecOfTrks[1]);
                        vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = myVertexAdaptive.trackWeight(vecOfTrks[2]);
                        vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = myVertexAdaptive.trackWeight(vecOfTrks[3]);
                    } else {
                        vtxcomp_adpt_isvalid[vtxcomp_n] = 0;
                        vtxcomp_adpt_x2dof[vtxcomp_n] = -9999;
                        vtxcomp_adpt_ndof[vtxcomp_n] = -9999;
                        vtxcomp_adpt_x2prob[vtxcomp_n] = -9999;
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][0]  = -9999;
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][1]  = -9999;
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][2]  = -9999;
                        vtxcomp_adpt_tk_chi2[vtxcomp_n][3]  = -9999;
                        vtxcomp_adpt_tk_weight[vtxcomp_n][0]  = -9999;
                        vtxcomp_adpt_tk_weight[vtxcomp_n][1]  = -9999;
                        vtxcomp_adpt_tk_weight[vtxcomp_n][2]  = -9999;
                        vtxcomp_adpt_tk_weight[vtxcomp_n][3]  = -9999;
                    }*/
                    vtxcomp_n++;

                } //end l
            } //end k
        } //end j
    } //end i



    return true;
}
