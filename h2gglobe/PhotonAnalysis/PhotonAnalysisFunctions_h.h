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
    TLorentzVector *calopos;
    bool pixSeed;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double sieie;
    double hoe;
};
