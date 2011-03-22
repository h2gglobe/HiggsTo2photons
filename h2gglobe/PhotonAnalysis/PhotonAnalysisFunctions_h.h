void PhotonAnalysisReducedOutputTree();
void myFillHistPhotonAnalysis(Util*, int );
void myFillHistPhotonAnalysisRed(Util * , int );
void myStatPhotonAnalysisRed(Util *, int);
void myReducePhotonAnalysis(Util * , int );
void myGetBranchPhotonAnalysis();
int myFillReducedVarPhotonAnalysis(Util * , int );
void mySetBranchAddressRedPhotonAnalysis();
int mySelectEventRedPhotonAnalysis(Util * , int);

void InitRealPhotonAnalysis(Util *, int typerun);
void TermRealPhotonAnalysis(int typerun);


//ADDED MJ
Float_t pho_Et[MAX_PHOTONS];

struct Elec{
    TLorentzVector *p4;
    TVector3 *calopos;
    bool pixSeed;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double sieie;
    double hoe;
};

//Added NCKW
bool ElecP4greater(Elec e1, Elec e2){
	
     return e1.p4->Pt() > e2.p4->Pt();


}
