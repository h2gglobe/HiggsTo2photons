// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::GlobeCtIsol(int mode, TLorentzVector* p4, float ptCut, float drCutMin, float drCutMax, Int_t & nIsol, Float_t & ptIsol, Float_t & angle1, Float_t & angle2, Float_t & angle3) {
  nIsol=0;
  ptIsol=0.;
  angle1=10.;
  angle2=10.;
  angle3=10.;

  //must put a track selection

  for (int i=0; i<ct_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) ct_p4->At(i);
    if(tempp4->Et()<ptCut) continue;
    double dr=p4->DeltaR(*tempp4);
    if(dr<drCutMin) continue;
    if(dr<angle1) {
      angle3=angle2;
      angle2=angle1;
      angle1=dr;
    }
    else if (dr<angle2) {
      angle3=angle2;
      angle2=dr;
    }
    else if (dr<angle3) {
      angle3=dr;
    }
    if(dr>drCutMax) continue;
    nIsol++;
    ptIsol+=tempp4->Et();
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::GlobeMatchIsl(TLorentzVector* p4, Float_t & deltaR) {
  deltaR=10.;
  int imatch=-1;

  for (int i=0; i<sc_islbar_n; i++) {
    TLorentzVector * tempp4= (TLorentzVector *) sc_islbar_p4->At(i);
    double dr=p4->DeltaR(*tempp4);
    if(dr<deltaR) {
      deltaR=dr;
      imatch=i;
    }
  }
  if(imatch==-1) {
    cout<<"ERROR GlobeMatchIsl found no match!!!! "<<endl;
  }
  else if(deltaR>0.3) {
    cout<<"Strange, GlobeMatchIsl deltaR="<<deltaR<<" etapho "<<p4->Eta()<<endl;
  }
  return imatch;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
#include "eIDCuts.h"
std::pair<bool, bool> LoopAll::ElectronId(int index, eIDLevel type) { 

  std::pair<bool, bool> isoIDResult(true, true);

  //TLorentzVector* p4 = (TLorentzVector*)sc_p4->At(el_std_scind[index]);
  TLorentzVector* p4 = (TLorentzVector*)el_std_sc->At(index);
  float eta = fabs(p4->Eta());
  int eb = 0, bin = 0;
  float see = 0;

  if (p4->Et() < 20.) 
    bin = 2;
  else if (p4->Et() > 30.)
    bin = 0;
  else
    bin =1;

  //#ifndef CMSSW3
  //if (eta < 1.479) {
  //see = sc_sieie[el_std_scind[index]];
  //eb = 0;
  //} else {
  //eb = 1; 
  //see = bc_sieie[sc_bcseedind[el_std_scind[index]]];
  //}
  //#else
  if (eta < 1.479) {
    see = el_std_sieiesc[index];
    eb = 0;
  } else {
    eb = 1; 
    see = el_std_sieie[index];
  }
  //#endif

  float eseedopincor = el_std_eseedopin[index] + el_std_fbrem[index];

  if(el_std_fbrem[index]<0) 
    eseedopincor = el_std_eseedopin[index]; 

  //#ifndef CMSSW3
  //float sip = sipCalculator(index);
  //#else
  float sip = fabs(el_std_ip_gsf[index]);
  //#endif

  int cat = ElectronClassification(index);
  
  float corr_tk_iso   = el_std_tkiso03[  index];
  float corr_ecal_iso = el_std_ecaliso04[index];
  float corr_hcal_iso = el_std_hcaliso04[index];  

  corr_tk_iso   = corr_tk_iso  *pow(40/p4->Et(), 2); 
  corr_ecal_iso = corr_ecal_iso*pow(40/p4->Et(), 2);
  corr_hcal_iso = corr_hcal_iso*pow(40/p4->Et(), 2);
  
  if ((corr_tk_iso > cutisotk[bin][type][cat]) ||
      (corr_ecal_iso > cutisoecal[bin][type][cat]) ||
      (corr_hcal_iso > cutisohcal[bin][type][cat]))

    isoIDResult.first = false;
  
  if (el_std_fbrem[index] < -2) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hoe[index]  > cuthoe[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (see > cutsee[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (fabs(el_std_dphiin[index]) > cutdphi[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;  
  }                            

  if (fabs(el_std_detain[index]) > cutdeta[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (eseedopincor < cuteopin[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }                      

  if (sip > cutip[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  if (el_std_hp_expin[index]  > cutmishits[bin][type][cat]) {
    isoIDResult.second = false;
    return isoIDResult;
  }

  return isoIDResult;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
int LoopAll::ElectronClassification(int index) {
  TLorentzVector* p4 = (TLorentzVector*) el_std_sc->At(index);

  int cat = -1;
  float eta = fabs(p4->Eta());

  if (eta < 1.479) {       // BARREL
    if(el_std_fbrem[index]<0.12)
      cat=1;
    else if (el_std_eopin[index] < 1.2 && el_std_eopin[index] > 0.9) 
      cat=0;
    else 
      cat=2;
  } else {                // ENDCAP
    if(el_std_fbrem[index]<0.2)
      cat=4;
    else if (el_std_eopin[index] < 1.22 && el_std_eopin[index] > 0.82) 
      cat=3;
    else 
      cat=5;
  }

  return cat;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
Float_t LoopAll::sipCalculator(int index) {

  Float_t ip = 0;

  if (el_std_tkind[index] != -1) {
    TLorentzVector* tk = (TLorentzVector*)tk_p4->At(el_std_tkind[index]);
    TVector3* my_tk_pos = (TVector3*)tk_vtx_pos->At(el_std_tkind[index]);

    // FIXME to handle the case of multiple vertices

    if (vtx_std_n != 0) {
      TVector3* my_vtx_pos = (TVector3*)vtx_std_xyz->At(0);

      // this is d0 "corrected" for the vertex...
      ip = fabs((-(my_tk_pos->X()-my_vtx_pos->X())*tk->Y()+(my_tk_pos->Y()-my_vtx_pos->Y()) * tk->X())/tk->Pt());
    } else {
      ip = fabs((-(my_tk_pos->X())*tk->Y()+(my_tk_pos->Y()) * tk->X())/tk->Pt());
    }

  }

  return ip;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::eIDInfo(Int_t index, Int_t& iso_result, Int_t& id_result, Int_t eIDMaxLevel) {

  iso_result = 0;
  id_result = 0;

  // FIXME add GetEntry functions
  for(Int_t i=0; i<eIDMaxLevel; ++i) {
    std::pair<bool, bool> result = ElectronId(index, LoopAll::eIDLevel(i));

    if (result.first) 
      iso_result = i;

    if (result.second)
      id_result = i;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------
// Vertex Analysis
// ---------------------------------------------------------------------------------------------------------------------------------------------
class GlobeVertexInfo : public VertexInfoAdapter
{
public:
	GlobeVertexInfo(LoopAll &);
	
	virtual int nvtx() const    { return lo_.vtx_std_n; };
	virtual int ntracks() const { return lo_.tk_n; };
	
	virtual bool hasVtxTracks() const { return true; }
	virtual const unsigned short * vtxTracks(int ii) const { return &(*lo_.vtx_std_tkind)[ii][0]; };
	virtual int vtxNTracks(int ii) const { return lo_.vtx_std_ntks[ii]; };
	virtual const float * vtxTkWeights(int ii) const { return &(*lo_.vtx_std_tkweight)[ii][0]; };

	virtual float tkpx(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Px(); };
	virtual float tkpy(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Py(); };
	virtual float tkpz(int ii) const { return ((TVector3*)lo_.tk_p4->At(ii))->Pz(); };
	
	virtual float tkPtErr(int ii) const { return lo_.tk_pterr[ii]; };
	virtual int   tkVtxId(int ii) const { return -1; };

	virtual float tkWeight(int ii, int jj) const { return (*lo_.vtx_std_tkweight)[jj][ii]; };
	
	virtual float vtxx(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->X(); };
	virtual float vtxy(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Y(); };
	virtual float vtxz(int ii) const { return ((TVector3*)lo_.vtx_std_xyz->At(ii))->Z(); };

	virtual float tkd0(int ii, int jj) const { return 0.; } // FIXME
	virtual float tkd0Err(int ii, int jj) const { return 0.; };  // FIXME

	virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
	virtual float tkdzErr(int ii, int jj) const { return 0.; };  // FIXME

	virtual bool tkIsHighPurity(int ii) const { return ( lo_.tk_quality[ii] & (1<<2) ) >> 2; };

	virtual ~GlobeVertexInfo();
	
private:
	LoopAll & lo_;
};


// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::GlobeVertexInfo(LoopAll & lo) : lo_(lo) {};

// ---------------------------------------------------------------------------------------------------------------------------------------------
GlobeVertexInfo::~GlobeVertexInfo() {};


// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::vertexAnalysis(HggVertexAnalyzer & vtxAna, int p1, int p2)
{
        GlobeVertexInfo vinfo(*this); 
	PhotonInfo
	  pho1(p1,*((TVector3*)pho_calopos->At(p1)),((TLorentzVector*)pho_p4->At(p1))->Energy()),
	  pho2(p2,*((TVector3*)pho_calopos->At(p2)),((TLorentzVector*)pho_p4->At(p2))->Energy());
	vtxAna.analyze(vinfo,pho1,pho2);

}

// ---------------------------------------------------------------------------------------------------------------------------------------------
std::vector<int> LoopAll::vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, int p1, int p2, std::vector<std::string> & vtxVarNames)
{
	assert( p1 == vtxAna.pho1() && p2 == vtxAna.pho2() );

	// preselect vertices : all vertices
        std::vector<int> preselAll;
        for(int i=0; i<vtx_std_n ; i++) {
          preselAll.push_back(i); 
        }

	// conversions infos
	PhotonInfo pho1(p1,
			*((TVector3*)pho_calopos->At(p1)),
			*((TVector3*) bs_xyz),
			*((TVector3*) pho_conv_vtx->At(p1)),
			((TLorentzVector*)pho_p4->At(p1))->Energy(),
			pho_isEB[p1],
			pho_conv_ntracks[p1],
			pho_conv_validvtx[p1],
			pho_conv_chi2_probability[p1] ,
			pho_conv_eoverp[p1]
			);
	

	PhotonInfo pho2(p2,
			*((TVector3*)pho_calopos->At(p2)),
			*((TVector3*) bs_xyz),
			*((TVector3*) pho_conv_vtx->At(p2)),
			((TLorentzVector*)pho_p4->At(p2))->Energy(),
			pho_isEB[p2],
			pho_conv_ntracks[p2],
			pho_conv_validvtx[p2],
			pho_conv_chi2_probability[p2] ,
			pho_conv_eoverp[p2]
			);


		
        float zconv = 0; 
        float dzconv = 0;
        std::vector<int> preselConv;


        if ( (pho_r9[p1] <0.93 || pho_r9[p2] <0.93) && (pho1.isAConversion() || pho2.isAConversion()) )  {
	  
          if (pho1.isAConversion()  && !pho2.isAConversion() ){
            zconv  = vtxAnaFromConv.vtxZ(pho1);
            dzconv = vtxAnaFromConv.vtxdZ(pho1);
          }
	  
          if (pho2.isAConversion() && !pho1.isAConversion()){
            zconv  = vtxAnaFromConv.vtxZ(pho2);
            dzconv = vtxAnaFromConv.vtxdZ(pho2);
          }
	  
          if ( pho1.isAConversion() && pho2.isAConversion()){
            float z1  = vtxAnaFromConv.vtxZ(pho1);
            float dz1 = vtxAnaFromConv.vtxdZ(pho1);
            
            float z2  = vtxAnaFromConv.vtxZ(pho2);
            float dz2 = vtxAnaFromConv.vtxdZ(pho2);
            
            zconv  = sqrt ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
            dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
          }
	  
	  // preselect vertices : only vertices in a window zconv +/- dzconv

	  for(int i=0; i < vtx_std_n; i++) {
	    TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(i);
	    if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
              preselConv.push_back(i); 
          }
	  
        }
	

	// preselection 
	if ( preselConv.size()==0 )
          vtxAna.preselection(preselAll);
        else 
          vtxAna.preselection(preselConv);

	

	std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
	cout << "\n\nRanks product" << endl;
	cout << "best vertex " << rankprod[0] << endl;
	for(int ii=0; ii<vtx_std_n; ++ii) {
		int vtxrank = find(rankprod.begin(), rankprod.end(), ii) - rankprod.begin();
		cout << "vertx " << ii << " rank " << vtxrank << " " << vtxAna.ptbal(ii) << " " << vtxAna.ptasym(ii) << " " << vtxAna.logsumpt2(ii) << endl;
	}
	
	return rankprod;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
TLorentzVector LoopAll::get_pho_p4(int ipho, int ivtx)
{
	PhotonInfo p(ipho, *((TVector3*)pho_calopos->At(ipho)),((TLorentzVector*)pho_p4->At(ipho))->Energy());
	TVector3 * vtx = (TVector3*) vtx_std_xyz->At(ivtx);
	return p.p4( vtx->X(), vtx->Y(), vtx->Z() );
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
void LoopAll::fillTrackIsolation(float tkIso_ptmin,float tkIso_outerCone,float tkIso_innerCone,float tkIso_etaStripHalfW,float tkIso_dzmax,float tkIso_dxymax)
{
	pho_trksumpthgg->clear(); pho_trksumpthgg->resize(pho_n,std::vector<float>(vtx_std_n,0.));
	pho_ntrkhgg->clear(); pho_ntrkhgg->resize(pho_n,std::vector<int>(vtx_std_n,0.));
	
	for(int ipho=0;ipho<pho_n;++ipho){
		for(int ivtx=0;ivtx<vtx_std_n;++ivtx){
			TLorentzVector p4 = get_pho_p4( ipho, ivtx );
			std::pair<Float_t,Int_t> tkIso = TrackIsoHgg(&p4, ivtx, tkIso_ptmin, tkIso_outerCone, tkIso_innerCone, 
								     tkIso_etaStripHalfW, tkIso_dzmax, tkIso_dxymax );
			(*pho_trksumpthgg)[ipho][ivtx] = tkIso.first;
			(*pho_ntrkhgg)[ipho][ivtx] = tkIso.second;
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------------------------------
std::pair<Float_t,Int_t> LoopAll::TrackIsoHgg(TLorentzVector *photon_p4, Int_t vtxind, Float_t PtMin, Float_t OuterConeRadius, 
					      Float_t InnerConeRadius, Float_t EtaStripHalfWidth, Float_t dzmax, Float_t dxymax) 
{
	// TRACKER Isolation
	if(vtxind<0)return std::make_pair(-99,-99);
	TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(vtxind);
	float SumTrackPt=0;
	int nTks = 0;
	for(unsigned int itk=0; itk!=tk_n; itk++) {
		TLorentzVector * tkp4= (TLorentzVector *) tk_p4->At(itk);
		if(tkp4->Pt() < PtMin)continue;
		TVector3 * tkpos= (TVector3 *) tk_vtx_pos->At(itk);
		double deltaz = fabs(vtxpos->Z() - tkpos->Z());
		if(deltaz > dzmax)continue;
		double dxy = ( -(tkpos->X() - vtxpos->X())*tkp4->Py() + (tkpos->Y() - vtxpos->Y())*tkp4->Px()) / tkp4->Pt();
		if(fabs(dxy) > dxymax)continue;

		double tk_eta = tkp4->Eta();
		double tk_phi = tkp4->Phi();
		double deta = fabs(photon_p4->Eta() - tk_eta);
		double dphi = fabs(photon_p4->Phi() - tk_phi);
		if(dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
		double deltaR = sqrt(deta*deta + dphi*dphi);
		if(deltaR < OuterConeRadius && deltaR >= InnerConeRadius && deta >= EtaStripHalfWidth) { 
			SumTrackPt+=tkp4->Pt();
			++nTks;
		}
	}

	return std::make_pair(SumTrackPt,nTks);
}
