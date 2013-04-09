#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeReducedGen.h"

GlobeReducedGen::GlobeReducedGen(const edm::ParameterSet& iConfig) {

  dR_min_for_matching = iConfig.getParameter<double>("GlobeReducedGendRMin");
}

void GlobeReducedGen::defineBranch(TTree* tree) {

  lptgeninfo_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
  lptgen_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
  lptgen_befrad_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
   
  tree->Branch("lptgeninfo_n", &lptgeninfo_n, "lptgeninfo_n/I");
  tree->Branch("lptgen_n", &lptgen_n, "lptgen_n/I");
  
  tree->Branch("lptgeninfo_p4", "TClonesArray", &lptgeninfo_p4, 32000, 0);
  tree->Branch("lptgen_p4", "TClonesArray", &lptgen_p4, 32000, 0);
  tree->Branch("lptgen_befrad_p4", "TClonesArray", &lptgen_befrad_p4, 32000, 0);
  
  tree->Branch("lptgeninfo_status", lptgeninfo_status, "lptgeninfo_status[lptgeninfo_n]/I");
  tree->Branch("lptgeninfo_pdgid", lptgeninfo_pdgid, "lptgeninfo_pdgid[lptgeninfo_n]/I");
  tree->Branch("lptgeninfo_mother", lptgeninfo_mother, "lptgeninfo_mother[lptgeninfo_n]/I");

  tree->Branch("lptgen_status", lptgen_status,           "lptgen_status[lptgen_n]/I"     );
  tree->Branch("lptgen_pdgid", lptgen_pdgid,             "lptgen_pdgid[lptgen_n]/I"      );
  tree->Branch("lptgen_mother", lptgen_mother,           "lptgen_mother[lptgen_n]/I"     );
  tree->Branch("lptgen_motherpdgid", lptgen_motherpdgid, "lptgen_motherpdgid[lptgen_n]/I");

  tree->Branch("lptgen_indrec", lptgen_indrec,           "lptgen_indrec[lptgen_n]/I");
  tree->Branch("lptgen_indrecel", lptgen_indrecel,       "lptgen_indrecel[lptgen_n]/I");
  tree->Branch("lptgen_indrecph", lptgen_indrecph,       "lptgen_indrecph[lptgen_n]/I");
  tree->Branch("lptgen_indrecmu", lptgen_indrecmu,       "lptgen_indrecmu[lptgen_n]/I");

  tree->Branch("lptgen_drmatch", lptgen_drmatch,         "lptgen_drmatch[lptgen_n]/F");
  tree->Branch("lptgen_drmatchel", lptgen_drmatchel,     "lptgen_drmatchel[lptgen_n]/F");
  tree->Branch("lptgen_drmatchmu", lptgen_drmatchmu,     "lptgen_drmatchmu[lptgen_n]/F");
  tree->Branch("lptgen_drmatchph", lptgen_drmatchph,     "lptgen_drmatchph[lptgen_n]/F");

  tree->Branch("lptgen_indinfo", lptgen_indinfo,         "lptgen_indinfo[lptgen_n]/I"    );
  tree->Branch("lptgen_historycode", lptgen_historycode, "lptgen_historycode[lptgen_n]/I");
}

void GlobeReducedGen::fillRedGenList(GlobeGenerator * gen, GlobeLeptons * lep){
  
  lptgeninfo_n=0;
  //lptgeninfo_p4->Clear();
  lptgen_n=0;
  // lptgen_p4->Clear();
  // lptgen_befrad_p4->Clear();
  
  int gen_keep[MAX_GENERATOR];
  int gen_historycode[MAX_GENERATOR];
  int temp_historycode[MAX_GENERATOR];
  
  //int newmother[MAX_GENERATOR];
  //int newline[MAX_GENERATOR];
  
  for (int i=0; i<gen->gen_n; i++) {
    gen_keep[i]=0;
    gen_historycode[i]=0;
    temp_historycode[i]=0;
  }
  
  for (int i=0; i<gen->gen_n; i++) {
    if(gen->gen_status[i]==1) {
      if(gen->gen_pdgid[i]==11
         ||gen->gen_pdgid[i]==-11
         ||gen->gen_pdgid[i]==13
         ||gen->gen_pdgid[i]==-13
         //||gen->gen_pdgid[i]==15 
         //||gen->gen_pdgid[i]==-15
         ||gen->gen_pdgid[i]==22) {
        
        int history[10];
        int nhistory;
        //TLorentzVector * pp4= (TLorentzVector *) gen->gen_p4->At(i);
        int rad_photons[10];
        int nrad_photons;
        int iselected = LeptonsGenInfo(gen, i, gen_keep, nhistory, history, nrad_photons, rad_photons);
        
        if(iselected)
          { gen_keep[i]=1; } // std::cout<<"I was selected"<<std::endl;}
        
        if(iselected) {
          //std::cout<<"nhistory "<<nhistory<<std::endl;
          for (int j=0; j<nhistory; j++) {
            //std::cout<<"history "<<j<<" "<<history[j]<<std::endl;
            
            //tag leptons
            if( fabs(history[j]) == 15 ) gen_historycode[i] = 15;
            if( fabs(history[j]) == 13 ) gen_historycode[i] = 13;
            if( fabs(history[j]) == 11 ) gen_historycode[i] = 11;
            
            //tag leptons from Zs
            if( fabs(history[j]) == 23 && gen_historycode[i] == 11) gen_historycode[i] = 1123;
            if( fabs(history[j]) == 23 && gen_historycode[i] == 13) gen_historycode[i] = 1323;
            if( fabs(history[j]) == 23 && gen_historycode[i] == 15) gen_historycode[i] = 1523;
            
            //tag leptons from Ws
            if( fabs(history[j]) == 24 && gen_historycode[i] == 11) gen_historycode[i] = 1124;
            if( fabs(history[j]) == 24 && gen_historycode[i] == 13) gen_historycode[i] = 1324;
            if( fabs(history[j]) == 24 && gen_historycode[i] == 15) gen_historycode[i] = 1524;
            
            //tag Zs from higgs
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1123) gen_historycode[i] = 2123;
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1323) gen_historycode[i] = 2323;
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1523) gen_historycode[i] = 2523;
            
            //tag Ws from higgs
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1124) gen_historycode[i] = 2124;
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1324) gen_historycode[i] = 2324;
            if( fabs(history[j]) == 25 && gen_historycode[i] == 1524) gen_historycode[i] = 2524;
            
          }
          //std::cout<<"nrad_photons "<<nrad_photons<<std::endl;
          //for (int j=0; j<nrad_photons; j++) {
          //std::cout<<"rad_photons "<<j<<" "<<rad_photons[j]<<std::endl;
          //}
        }
      }
    }
  }
  
  
  //FILL REDUCED INFO
  
  int gen_newid[MAX_GENERATOR];
  
  //#define MAX_LPT_GENINFO 100
  //lptgeninfo_n=0;
  for (int i=0; i<gen->gen_n; i++) {
    if(lptgeninfo_n<MAX_LPT_GENINFO) {
      lptgeninfo_mother[lptgeninfo_n]=-1;
      if(gen_keep[i]) {
        gen_newid[i]=lptgeninfo_n;
        TLorentzVector * pp4= (TLorentzVector *) gen->gen_p4->At(i);
        new ((*lptgeninfo_p4)[lptgeninfo_n]) TLorentzVector(*pp4);
        lptgeninfo_status[lptgeninfo_n]=gen->gen_status[i];
        lptgeninfo_pdgid[lptgeninfo_n]=gen->gen_pdgid[i];
        temp_historycode[lptgeninfo_n]=gen_historycode[i];
        
        if(gen->gen_mother[i]>=0) {
          lptgeninfo_mother[lptgeninfo_n]=gen_newid[gen->gen_mother[i]];
        }
        else {
          lptgeninfo_mother[lptgeninfo_n]=-1;
        }
        lptgeninfo_n++;
      }
      else {
        gen_newid[i]=-1;
      }
    } else {
      std::cout<<"ATTENTION too many MAX_LPT_GENINFO"<<std::endl;
      break;
    }
  }
  
  //FINISHED GETTING GEN INFO
  
  //cuts:
  //CHECK THESE ARE HARDCODED
  float leptons_cut_et_gen=0.0;
  float leptons_cut_eta_gen=3.3;
  
  //MAKE LEPTON GEN LIST FROM LPTGENINFO
  for (int i=0; i<lptgeninfo_n; i++) {
    if(lptgeninfo_status[i]==1) {
      TLorentzVector * pp4= (TLorentzVector *) lptgeninfo_p4->At(i);
      if(pp4->Pt()<leptons_cut_et_gen) continue;
      
      if(fabs((float) pp4->Eta())>leptons_cut_eta_gen) continue;
      
      if(lptgen_n>=MAX_LPT_GEN) continue;
      new ((*lptgen_p4)[lptgen_n]) TLorentzVector(*pp4);
      lptgen_pdgid[lptgen_n]=lptgeninfo_pdgid[i];
      lptgen_indinfo[lptgen_n]=i;
      
      //this was wrong 
      //lptgen_motherpdgid[lptgen_n]=gen->gen_pdgid[gen->gen_mother[i]];//-100; //not filled
      
      lptgen_mother[lptgen_n] = lptgeninfo_mother[i];
      lptgen_motherpdgid[lptgen_n] =
        (lptgeninfo_mother[i] < 0 ? -100 : lptgeninfo_pdgid[lptgeninfo_mother[i]]);// -100 was not filled
      
      lptgen_indrec[lptgen_n]=-1;
      lptgen_indrecel[lptgen_n]=-1;
      lptgen_indrecmu[lptgen_n]=-1;
      lptgen_indrecph[lptgen_n]=-1;
      lptgen_status[lptgen_n]=lptgeninfo_status[i];
      lptgen_drmatch[lptgen_n]=10.;
      lptgen_drmatchel[lptgen_n]=10.;
      lptgen_drmatchmu[lptgen_n]=10.;
      lptgen_drmatchph[lptgen_n]=10.;
           
      lptgen_historycode[lptgen_n]=temp_historycode[i];//0; //not filled, to be made a function

      //even a tau may radiate a photon
      if(
         (lptgeninfo_pdgid[i]==11 ||
          lptgeninfo_pdgid[i]==-11 ||
          lptgeninfo_pdgid[i]==13 ||
          lptgeninfo_pdgid[i]==-13) && lptgeninfo_mother[i]>=0) {
        
        if(lptgeninfo_pdgid[lptgeninfo_mother[i]]==lptgeninfo_pdgid[i]) {
          TLorentzVector * pp4_new= (TLorentzVector *) lptgeninfo_p4->At(lptgeninfo_mother[i]);
          new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4_new);
        }
        else {
          new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4);
        }
      }
      else {
        new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4);
      }
      lptgen_n++;
    }
  }
  //now do the matching gen rec between these two
  //missing lptgen_indrec set to -1;
  
  for (int i=0; i<lep->lpt_n; i++) {
    TLorentzVector * pp4= (TLorentzVector *) lep->lpt_p4->At(i);

    float deltaR=dR_min_for_matching;

    //imatch = GlobeMatchWithGen(gen, pp4, deltaR, lep->lpt_pdgid[i], 3, 0.2);
    for(int type = 0; type<3; type++) {
      int imatch = -1;
      
      if (type == 0) {
        imatch = GlobeMatchWithGen(gen, pp4, deltaR, 11, 3, 0.2);
        if(imatch!=-1) { 
          if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
          }
          //added this for cross matching, should also add the same variables for the reverse matching
          else if (abs(lep->lpt_pdgid[i]) == 13) {
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
          }
          else if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
          }
        }
      }
      
      if (type == 1) {
        imatch = GlobeMatchWithGen(gen, pp4, deltaR, 13, 3, 0.2);
        if(imatch!=-1) { 
          if (abs(lep->lpt_pdgid[i]) == 13) {
            lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
          }     
          else if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
          }
          else if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
          }
        }
      }
      
      if (type == 2) {
        imatch = GlobeMatchWithGen(gen, pp4, deltaR, 22, 3, 0.2);
        if(imatch!=-1) { 
          if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
          }
          else if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
          }
          else if (abs(lep->lpt_pdgid[i]) == 13) {
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
          }
        }
      }      
    }
  }
}

int GlobeReducedGen::LeptonsGenInfo(GlobeGenerator * gen, int j, int * gen_keep, int & nhistory, int * history, int & nrad_photons, int * rad_photons) {
  //, int & imothpdg, int & iquarkmoth, int & idquarkmoth, int & iwzmoth, int & idwzmoth, int & ihiggsmoth, int & idhiggsmoth) {
  nhistory=0;
  nrad_photons=0;
  //int history[10];
  
  int pdgidMother=-1;
  
  //int j=i;
  int motherId=-1;
  int pdgidPrevious=-1;
  
  int printit=0;
  
  int nloop=0;
  
  
  //MARCO TRIED TO FIX IT, I AM NOT SURE
  //while(j!=-1&&j!=0) {//attention to j!=0 
  while(j!=-1) {//attention to j!=0 
    
    //std::cout<<"here 1 "<<j<<std::endl;
    
    pdgidPrevious=gen->gen_pdgid[j];
    if(pdgidPrevious<0) pdgidPrevious=-pdgidPrevious;
    
    motherId=gen->gen_mother[j];
    
    //std::cout<<"here motherId "<<motherId<<std::endl;
    
    if(motherId>=0) {
      
      pdgidMother=0;
      int pdgid=gen->gen_pdgid[motherId];
      if(pdgid<0) pdgid=-pdgid;
      int pdgidfull=pdgid;
      if(pdgid>10000) pdgid=pdgid%10000;

      if(pdgidPrevious==22 && (pdgid==21 || (pdgid>-26 && pdgid<26))) {
        pdgidMother=pdgid;
      }
      else if(pdgidPrevious!=22 && pdgid>1000) {
        pdgidMother=pdgid/1000;
      }
      else if(pdgidPrevious!=22 && pdgid>100) {
        pdgidMother=pdgid/100;
      }
      else if (pdgidPrevious!=22 && (pdgid==15 || pdgid==13 || pdgid==11 || (pdgid>=22 && pdgid<=25) || pdgid==6)) { 
        pdgidMother=pdgid;
      }
      else {
        if(pdgidPrevious!=22 && pdgidfull>6 && pdgidfull!=21 && pdgidfull!=92 && pdgidfull!=91) { //ignore the quarks for now
          //std::cout<<"LeptonsGenInfo: Attention HERE TO BE CONSIDERED pdgidPrevious "<<pdgidPrevious<<" pdgidfull "<<pdgidfull<<" motherId "<<motherId<<" gen_status[motherId] "<<gen->gen_status[motherId]<<" gen_pdgid[motherId] "<<gen->gen_pdgid[motherId]<<" pdgidMother "<<pdgidMother<<std::endl;
        }
      }

      if(pdgidMother<4) break;

      gen_keep[motherId]=1;

      int in=0;
      if(pdgidMother>0) {
        if(pdgidMother==pdgidPrevious) { //this should be only leptons with radiative photons
          int imax=j+30;
          if(gen->gen_n>imax) imax=gen->gen_n;
          for(int k=motherId; k<imax; k++) {
            if(gen->gen_mother[k]==motherId && k!=j ) {
              if(gen->gen_pdgid[k]==22) {

                //THESE ARE THE PHOTONS TO BE KEPT SOMEWHERE
                //TLorentzVector * pp4= (TLorentzVector *) gen->gen_p4->At(k);
                //if(pp4->Pt()>0.1) 
                {
                  if(in>-1) { //just for the compiler not to complain
                    //std::cout<<"DAUG only photons pdgidMother "<<pdgidPrevious<<" in "<<in++<<" "<<gen->gen_status[k]<<" pdgid daugh "<<gen->gen_pdgid[k]<<" "<<pp4->Pt()<<std::endl;
                    //keep the photons
                    gen_keep[k]=1;
                    rad_photons[nrad_photons++]=k;  
                  }
                }
              }
            }
          }
        }

        if(in>0) {
          //printit=1;
        }

        //add the particle in the history
        if(nhistory<10) {
          if(nhistory==0) 
            history[nhistory++]=pdgidPrevious;
          if(pdgidMother!=0) {
            if(pdgidMother!=history[nhistory-1] || pdgidMother==11 || pdgidMother==13 || pdgidMother==15 ) 
              //if(pdgidMother!=history[nhistory-1]) 
              {
                history[nhistory++]=pdgidMother;
              }
            //if(pdgidMother==history[nhistory-1])  //4,5,6,23,25 all the time
            //std::std::cout<<"mother like daughter: "<<pdgidMother<<std::endl;
          }
        }
      }
    }
    else {
      //here it should be kept for single particles (MARCO)
      //std::cout<<"here nloop "<<nloop<<std::endl;
      if(nloop==0) {
	history[nhistory++]=pdgidPrevious;
	//std::cout<<"here nhistory "<<nhistory<<std::endl;
	return 1;
      }
    }
    j=gen->gen_mother[j];
    nloop++;
  }

  //at the end prints the history
  if(printit) {
    std::cout<<"nhistory "<<nhistory<<std::endl;
    for (int j=0; j<nhistory; j++) {
      std::cout<<"history "<<j<<" "<<history[j]<<std::endl;
    }
  }

  return nhistory;
}



int  GlobeReducedGen::GlobeMatchWithGen(GlobeGenerator * gen, TLorentzVector* p4, Float_t & deltaR, int pdgid, int gencoll, float cutgenrecdrmatch)
{
  
  //gencoll 1 full, 2 reduced leptons, 3 reduced list of leptons (FOR NOW USE 2 or 3)
  
  deltaR=10.;
  //float deltaPhi=10.;
  //float deltaEta=10.;
  int imatch=-1;
  
  int n_n=gen->gen_n;
  if(gencoll==2) n_n=lptgeninfo_n;
  if(gencoll==3) n_n=lptgen_n;
  //std::cout<<"here"<<std::endl;
  
  for (int i=0; i<n_n; i++) {
    //attention to the right photon matching!!!
    //int j=-10000000;
    //std::cout<<"here"<<i<<" "<<lptgen_status[i]<<" "<<lptgen_pdgid[i]<<std::endl;
    
    if(gencoll==1) {
      if(gen->gen_status[i]!=1) continue;
      if(gen->gen_pdgid[i]!=pdgid && gen->gen_pdgid[i]!=-pdgid) continue;
      
      //if(gen->gen_mother[i]>=0) {
      //  j=gen->gen_pdgid[gen->gen_mother[i]];
      //}
      //ignore those with wrong mother NOT IMPLEMENTED YET
      //if(j!=21 && j!=22 && (j<-20 || j>20)) continue;
    }
    else if(gencoll==2) {
      if(lptgeninfo_status[i]!=1) continue;
      if(lptgeninfo_pdgid[i]!=pdgid && lptgeninfo_pdgid[i]!=-pdgid) continue;
    }
    else if(gencoll==3) {
      if(lptgen_status[i]!=1) continue;
      if(lptgen_pdgid[i]!=pdgid && lptgen_pdgid[i]!=-pdgid) continue;
    }
    
    //std::cout<<"here"<<std::endl;
    
    TLorentzVector * tempp4= (TLorentzVector *) gen->gen_p4->At(i);
    if(gencoll==2) {
      tempp4= (TLorentzVector *) lptgeninfo_p4->At(i);
    }
    if(gencoll==3) {
      tempp4= (TLorentzVector *) lptgen_p4->At(i);
    }
    double dr=p4->DeltaR(*tempp4);
    //double dphi=p4->DeltaPhi(*tempp4);
    //float deta=fabs((float) p4->Eta()- tempp4->Eta());
    if(dr<deltaR) {
      deltaR=dr;
      //deltaPhi=dphi;
      //deltaEta=deta;
      imatch=i;
    }
  }
  //if(deltaR>0.01 && deltaR<cutgenrecdrmatch) std::cout<<"dr, dphi, deta "<<deltaR<<" "<<deltaPhi<<" "<<deltaEta<<std::endl;
  if(deltaR>cutgenrecdrmatch)
    imatch=-1;
  return imatch;
}
/*
  int LoopAll::GlobeGenHistoryCode(int i, int gencoll) {
  //gencoll = 0 gen, 1, lptgeninfo, 2 lptgen
  int ipart=-1;
  if(gencoll==2) {
  ipart=lptgen_indinfo[i];
  }
  return 0;
  //to be finished
  
  }
*/



void GlobeReducedGen::fillRedGenList(GlobeGenParticles* gp, GlobeLeptons * lep){
  
  lptgeninfo_n=0;
  //lptgeninfo_p4->Clear();
  lptgen_n=0;
  // lptgen_p4->Clear();
  // lptgen_befrad_p4->Clear();
  
  int gp_keep[MAX_GENERATOR];
  int gp_historycode[MAX_GENERATOR];
  int temp_historycode[MAX_GENERATOR];
  
  //int newmother[MAX_GENERATOR];
  //int newline[MAX_GENERATOR];
  
  for (int i=0; i<gp->gp_n; i++) {
    gp_keep[i]=0;
    gp_historycode[i]=0;
    temp_historycode[i]=0;
  }
  
  for (int i=0; i<gp->gp_n; i++) {

    if(gp->gp_status[i] == 1) {

      if(gp->gp_pdgid[i] == 11 or gp->gp_pdgid[i] == -11 or
         gp->gp_pdgid[i] == 13 or gp->gp_pdgid[i] == -13 or
         gp->gp_pdgid[i] == 22) {
        
        int history[10];
        int nhistory;
        //TLorentzVector * pp4= (TLorentzVector *) gp->gp_p4->At(i);
        int rad_photons[10];
        int nrad_photons;
        if (LeptonsGenInfo(gp, i, gp_keep, nhistory, history, nrad_photons, rad_photons)) {
	  gp_keep[i] = 1;

          //std::cout<<"nhistory "<<nhistory<<std::endl;
          for (int j=0; j<nhistory; j++) {
            //std::cout<<"history "<<j<<" "<<history[j]<<std::endl;

            //tag leptons
            if( fabs(history[j]) == 15 ) gp_historycode[i] = 15;
            if( fabs(history[j]) == 13 ) gp_historycode[i] = 13;
            if( fabs(history[j]) == 11 ) gp_historycode[i] = 11;

            //tag leptons from Zs
            if( fabs(history[j]) == 23 && gp_historycode[i] == 11) gp_historycode[i] = 1123;
            if( fabs(history[j]) == 23 && gp_historycode[i] == 13) gp_historycode[i] = 1323;
            if( fabs(history[j]) == 23 && gp_historycode[i] == 15) gp_historycode[i] = 1523;

            //tag leptons from Ws
            if( fabs(history[j]) == 24 && gp_historycode[i] == 11) gp_historycode[i] = 1124;
            if( fabs(history[j]) == 24 && gp_historycode[i] == 13) gp_historycode[i] = 1324;
            if( fabs(history[j]) == 24 && gp_historycode[i] == 15) gp_historycode[i] = 1524;

            //tag Zs from higgs
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1123) gp_historycode[i] = 2123;
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1323) gp_historycode[i] = 2323;
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1523) gp_historycode[i] = 2523;

            //tag Ws from higgs
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1124) gp_historycode[i] = 2124;
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1324) gp_historycode[i] = 2324;
            if( fabs(history[j]) == 25 && gp_historycode[i] == 1524) gp_historycode[i] = 2524;
            
          }
        }
      }
    }
  }
  











  
  //FILL REDUCED INFO
  
  int gp_newid[MAX_GENERATOR];
  
  //#define MAX_LPT_GENINFO 100
  //lptgeninfo_n=0;
  std::cout << lptgeninfo_n << " " << gp->gp_n << std::endl;

  for (int i=0; i<gp->gp_n; i++) {
    if(lptgeninfo_n<MAX_LPT_GENINFO) {
      lptgeninfo_mother[lptgeninfo_n]=-1;
      if(gp_keep[i]) {
        gp_newid[i]=lptgeninfo_n;
        TLorentzVector * pp4= (TLorentzVector *) gp->gp_p4->At(i);
        new ((*lptgeninfo_p4)[lptgeninfo_n]) TLorentzVector(*pp4);
        lptgeninfo_status[lptgeninfo_n]=gp->gp_status[i];
        lptgeninfo_pdgid[lptgeninfo_n]=gp->gp_pdgid[i];
        temp_historycode[lptgeninfo_n]=gp_historycode[i];
	 
        if(gp->gp_mother[i]>=0) {
          lptgeninfo_mother[lptgeninfo_n]=gp_newid[gp->gp_mother[i]];
        }
        else {
          lptgeninfo_mother[lptgeninfo_n]=-1;
        }
        lptgeninfo_n++;
      }
      else {
        gp_newid[i]=-1;
      }
    }
    else {
      std::cout<<"ATTENTION too many MAX_LPT_GENINFO"<<std::endl;
      break;
    }
  }
  
  //FINISHED GETTING GEN INFO
  
  //cuts:
  //CHECK THESE ARE HARDCODED
  float leptons_cut_et_gen=0.0;
  float leptons_cut_eta_gen=3.3;

  //MAKE LEPTON GEN LIST FROM LPTGENINFO
  for (int i=0; i<lptgeninfo_n; i++) {
    if(lptgeninfo_status[i]==1) {
      TLorentzVector * pp4= (TLorentzVector *) lptgeninfo_p4->At(i);
      if(pp4->Pt()<leptons_cut_et_gen) continue;

      if(fabs((float) pp4->Eta())>leptons_cut_eta_gen) continue;
     
      if(lptgen_n>=MAX_LPT_GEN) continue;
      new ((*lptgen_p4)[lptgen_n]) TLorentzVector(*pp4);
      lptgen_pdgid[lptgen_n]=lptgeninfo_pdgid[i];
      lptgen_indinfo[lptgen_n]=i;

      //this was wrong 
      //lptgen_motherpdgid[lptgen_n]=gen->gen_pdgid[gen->gen_mother[i]];//-100; //not filled

      if(lptgeninfo_mother[i]<0) {
	lptgen_motherpdgid[lptgen_n]=-100;
	lptgen_mother[lptgen_n] = lptgeninfo_mother[i];
      }
      else { 
	if(lptgeninfo_pdgid[lptgeninfo_mother[i]]!= lptgeninfo_pdgid[i] || lptgeninfo_status[lptgeninfo_mother[i]]!=3) {	
	  lptgen_motherpdgid[lptgen_n] = lptgeninfo_pdgid[lptgeninfo_mother[i]];
	  lptgen_mother[lptgen_n] = lptgeninfo_mother[i];
	}
	else {
	  if(lptgeninfo_mother[lptgeninfo_mother[i]]<0) {
	    lptgen_motherpdgid[lptgen_n]=-100;
	    lptgen_mother[lptgen_n]=lptgeninfo_mother[lptgeninfo_mother[i]];
	  }
	  else {
	    lptgen_motherpdgid[lptgen_n]=lptgeninfo_pdgid[lptgeninfo_mother[lptgeninfo_mother[i]]];
	    lptgen_mother[lptgen_n]=lptgeninfo_mother[lptgeninfo_mother[i]];
	  }
	}
      }

      /* OLD WRONG FOR ESTAR
      lptgen_mother[lptgen_n] = lptgeninfo_mother[i];
      lptgen_motherpdgid[lptgen_n] =
	(lptgeninfo_mother[i] < 0 ?
	 -100 :
	 lptgeninfo_pdgid[lptgeninfo_mother[i]]);// -100 was not filled
      */


      lptgen_indrec[lptgen_n]=-1;
      lptgen_indrecel[lptgen_n]=-1;
      lptgen_indrecmu[lptgen_n]=-1;
      lptgen_indrecph[lptgen_n]=-1;
      lptgen_status[lptgen_n]=lptgeninfo_status[i];
      lptgen_drmatch[lptgen_n]=10.;
      lptgen_drmatchel[lptgen_n]=10.;
      lptgen_drmatchmu[lptgen_n]=10.;
      lptgen_drmatchph[lptgen_n]=10.;
      
           
      lptgen_historycode[lptgen_n]=temp_historycode[i];//0; //not filled, to be made a function

      //even a tau may radiate a photon
      if(
         (lptgeninfo_pdgid[i]==11 ||
          lptgeninfo_pdgid[i]==-11 ||
          lptgeninfo_pdgid[i]==13 ||
          lptgeninfo_pdgid[i]==-13) && lptgeninfo_mother[i]>=0) {

        if(lptgeninfo_pdgid[lptgeninfo_mother[i]]==lptgeninfo_pdgid[i]) {
          TLorentzVector * pp4_new= (TLorentzVector *) lptgeninfo_p4->At(lptgeninfo_mother[i]);
          new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4_new);
        }
        else {
          new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4);
        }
      }
      else {
        new ((*lptgen_befrad_p4)[lptgen_n]) TLorentzVector(*pp4);
      }
      lptgen_n++;
    }
  }
  //now do the matching gen rec between these two
  //missing lptgen_indrec set to -1;
 
  for (int i=0; i<lep->lpt_n; i++) {
    TLorentzVector * pp4= (TLorentzVector *) lep->lpt_p4->At(i);

    float deltaR=dR_min_for_matching;

    //imatch = GlobeMatchWithGen(gen, pp4, deltaR, lep->lpt_pdgid[i], 3, 0.2);
    for(int type = 0; type<3; type++) {
      int imatch = -1;
      
      if (type == 0) {
        imatch = GlobeMatchWithGen(gp, pp4, deltaR, 11, 3, 0.2);
	if(imatch!=-1) { 
	  if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
          }
	  //added this for cross matching, should also add the same variables for the reverse matching
	  else if (abs(lep->lpt_pdgid[i]) == 13) {
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
	  }
	  else if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
	  }
        }
      }

      if (type == 1) {
        imatch = GlobeMatchWithGen(gp, pp4, deltaR, 13, 3, 0.2);
	if(imatch!=-1) { 
	  if (abs(lep->lpt_pdgid[i]) == 13) {
	    lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
          }     
	  else if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
	  }
	  else if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
	  }
        }
      }
       
      if (type == 2) {
        imatch = GlobeMatchWithGen(gp, pp4, deltaR, 22, 3, 0.2);
	if(imatch!=-1) { 
	  if (abs(lep->lpt_pdgid[i]) == 22) {
            lptgen_drmatch[imatch]=deltaR;
            lptgen_indrec[imatch]=i;
            lep->lpt_indgen[i]=imatch;
            lep->lpt_drmatch[i]=deltaR;
            lptgen_drmatchph[imatch]=deltaR;
            lptgen_indrecph[imatch]=i;
          }
	  else if (abs(lep->lpt_pdgid[i]) == 11) {
            lptgen_drmatchel[imatch]=deltaR;
            lptgen_indrecel[imatch]=i;
	  }
	  else if (abs(lep->lpt_pdgid[i]) == 13) {
            lptgen_drmatchmu[imatch]=deltaR;
            lptgen_indrecmu[imatch]=i;
	  }
        }
      }

    }
  }
}

int GlobeReducedGen::LeptonsGenInfo(GlobeGenParticles* gp, int j, int * gp_keep, int & nhistory, int * history, int & nrad_photons, int * rad_photons) {

  nhistory=0;
  nrad_photons=0;
  int pdgidMother=-1;
  
  int motherId=-1;
  int pdgidPrevious=-1;
  int printit=0;
  int nloop=0;

  while(j != -1) {//attention to j!=0 
    pdgidPrevious = gp->gp_pdgid[j];
    if (pdgidPrevious < 0) 
      pdgidPrevious = -pdgidPrevious;

    motherId = gp->gp_mother[j];

    if (motherId >= 0) {
      pdgidMother = 0;
      int pdgid= gp->gp_pdgid[motherId];
      if(pdgid < 0) 
	pdgid=-pdgid;
      int pdgidfull = pdgid;
      if(pdgid > 10000) 
	pdgid=pdgid%10000;

      if(pdgidPrevious==22 && (pdgid==21 || (pdgid>-26 && pdgid<26))) {
        pdgidMother=pdgid;
      }
      else if(pdgidPrevious!=22 && pdgid>1000) {
        pdgidMother=pdgid/1000;
      }
      else if(pdgidPrevious!=22 && pdgid>100) {
        pdgidMother=pdgid/100;
      }
      else if (pdgidPrevious!=22 && (pdgid==15 || pdgid==13 || pdgid==11 || (pdgid>=22 && pdgid<=25) || pdgid==6)) { 
        pdgidMother=pdgid;
      }
      else {
        if(pdgidPrevious!=22 && pdgidfull>6 && pdgidfull!=21 && pdgidfull!=92 && pdgidfull!=91) { //ignore the quarks for now
          //std::cout<<"LeptonsGenInfo: Attention HERE TO BE CONSIDERED pdgidPrevious "<<pdgidPrevious<<" pdgidfull "<<pdgidfull<<" motherId "<<motherId<<" gen_status[motherId] "<<gen->gen_status[motherId]<<" gen_pdgid[motherId] "<<gen->gen_pdgid[motherId]<<" pdgidMother "<<pdgidMother<<std::endl;
        }
      }

      if(pdgidMother<4) break;

      gp_keep[motherId]=1;

      int in=0;
      if(pdgidMother>0) {
        if(pdgidMother==pdgidPrevious) { //this should be only leptons with radiative photons
          int imax=j+30;
          if(gp->gp_n>imax) imax=gp->gp_n;
          for(int k=motherId; k<imax; k++) {
            if(gp->gp_mother[k]==motherId && k!=j ) {
              if(gp->gp_pdgid[k]==22) {

                //THESE ARE THE PHOTONS TO BE KEPT SOMEWHERE
                //TLorentzVector * pp4= (TLorentzVector *) gen->gen_p4->At(k);
                //if(pp4->Pt()>0.1) 
                {
                  if(in>-1) { //just for the compiler not to complain
                    //std::cout<<"DAUG only photons pdgidMother "<<pdgidPrevious<<" in "<<in++<<" "<<gen->gen_status[k]<<" pdgid daugh "<<gen->gen_pdgid[k]<<" "<<pp4->Pt()<<std::endl;
                    //keep the photons
                    gp_keep[k]=1;
                    rad_photons[nrad_photons++]=k;  
                  }
                }
              }
            }
          }
        }

        if(in>0) {
          //printit=1;
        }

        //add the particle in the history
        if(nhistory<10) {
          if(nhistory==0) 
            history[nhistory++]=pdgidPrevious;
          if(pdgidMother!=0) {
            if(pdgidMother!=history[nhistory-1] || pdgidMother==11 || pdgidMother==13 || pdgidMother==15 ) 
              //if(pdgidMother!=history[nhistory-1]) 
              {
                history[nhistory++]=pdgidMother;
              }
            //if(pdgidMother==history[nhistory-1])  //4,5,6,23,25 all the time
            //std::std::cout<<"mother like daughter: "<<pdgidMother<<std::endl;
          }
        }
      }
    }
    else {
      //here it should be kept for single particles (MARCO)
      //std::cout<<"here nloop "<<nloop<<std::endl;
      if(nloop==0) {
	history[nhistory++]=pdgidPrevious;
	//std::cout<<"here nhistory "<<nhistory<<std::endl;
	return 1;
      }
    }
    j=gp->gp_mother[j];
    nloop++;
  }

  //at the end prints the history
  if(printit) {
    std::cout<<"nhistory "<<nhistory<<std::endl;
    for (int j=0; j<nhistory; j++) {
      std::cout<<"history "<<j<<" "<<history[j]<<std::endl;
    }
  }

  return nhistory;
}



int  GlobeReducedGen::GlobeMatchWithGen(GlobeGenParticles* gp, TLorentzVector* p4, Float_t & deltaR, int pdgid, int gencoll, float cutgenrecdrmatch)
{
  
  //gencoll 1 full, 2 reduced leptons, 3 reduced list of leptons (FOR NOW USE 2 or 3)
  
  deltaR=10.;
  //float deltaPhi=10.;
  //float deltaEta=10.;
  int imatch=-1;
  
  int n_n=gp->gp_n;
  if(gencoll==2) n_n=lptgeninfo_n;
  if(gencoll==3) n_n=lptgen_n;
  //std::cout<<"here"<<std::endl;
  
  for (int i=0; i<n_n; i++) {
    //attention to the right photon matching!!!
    //int j=-10000000;
    //std::cout<<"here"<<i<<" "<<lptgen_status[i]<<" "<<lptgen_pdgid[i]<<std::endl;
    
    if(gencoll==1) {
      if(gp->gp_status[i]!=1) continue;
      if(gp->gp_pdgid[i]!=pdgid && gp->gp_pdgid[i]!=-pdgid) continue;
      
      //if(gp->gp_mother[i]>=0) {
      //  j=gp->gp_pdgid[gp->gp_mother[i]];
      //}
      //ignore those with wrong mother NOT IMPLEMENTED YET
      //if(j!=21 && j!=22 && (j<-20 || j>20)) continue;
    }
    else if(gencoll==2) {
      if(lptgeninfo_status[i]!=1) continue;
      if(lptgeninfo_pdgid[i]!=pdgid && lptgeninfo_pdgid[i]!=-pdgid) continue;
    }
    else if(gencoll==3) {
      if(lptgen_status[i]!=1) continue;
      if(lptgen_pdgid[i]!=pdgid && lptgen_pdgid[i]!=-pdgid) continue;
    }
    
    //std::cout<<"here"<<std::endl;
    
    TLorentzVector * tempp4= (TLorentzVector *) gp->gp_p4->At(i);
    if(gencoll==2) {
      tempp4= (TLorentzVector *) lptgeninfo_p4->At(i);
    }
    if(gencoll==3) {
      tempp4= (TLorentzVector *) lptgen_p4->At(i);
    }
    double dr=p4->DeltaR(*tempp4);
    //double dphi=p4->DeltaPhi(*tempp4);
    //float deta=fabs((float) p4->Eta()- tempp4->Eta());
    if(dr<deltaR) {
      deltaR=dr;
      //deltaPhi=dphi;
      //deltaEta=deta;
      imatch=i;
    }
  }
  //if(deltaR>0.01 && deltaR<cutgenrecdrmatch) std::cout<<"dr, dphi, deta "<<deltaR<<" "<<deltaPhi<<" "<<deltaEta<<std::endl;
  if(deltaR>cutgenrecdrmatch)
    imatch=-1;
  return imatch;
}
/*
  int LoopAll::GlobeGenHistoryCode(int i, int gencoll) {
  //gencoll = 0 gen, 1, lptgeninfo, 2 lptgen
  int ipart=-1;
  if(gencoll==2) {
  ipart=lptgen_indinfo[i];
  }
  return 0;
  //to be finished
  
  }
*/
