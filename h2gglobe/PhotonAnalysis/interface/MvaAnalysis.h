#ifndef __MVAANALYSIS__
#define __MVAANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include "MassResolution.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class MvaAnalysis : public PhotonAnalysis 
{
public:
	
	MvaAnalysis();
	virtual ~MvaAnalysis();
	
	virtual const std::string & name() const { return name_; };
	
	// LoopAll analysis interface implementation
	void Init(LoopAll&);
	void Term(LoopAll&);
	
	void GetBranches(TTree *, std::set<TBranch *>& );
	
	virtual bool SelectEvents(LoopAll&, int);
	virtual void ResetAnalysis();
	virtual void Analysis(LoopAll&, Int_t);

//	double GetPhotonResolution(LoopAll &,int);
//	double GetAngleResolutionCorrVtx(LoopAll &,int, int, int);
//	double GetAngleResolutionWrongVtx(LoopAll &, int, int, int);
//	double PropagateDz(LoopAll &, int, int, int, double);

//	double SecH(double);
//	double TanH(double);
	int SignalType(int);
	void SetBDTInputVariables(TLorentzVector*, TLorentzVector*, double, double, double, double, double, int cat = 0);
	void SetBDTInputTree(TTree *);
	
	// Options
	float leadEtCut;
	float subleadEtCut;
	std::string efficiencyFile;
	
	// EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	bool  doMCSmearing;
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
	bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear;
	float systRange;
	int   nSystSteps;   
	int   nEtaCategories, nR9Categories, nPtCategories;
	float massMin, massMax;
	int nDataBins;	
        float signalRegionWidth;
        float sidebandWidth;
        float sidebandShift;
        int numberOfSidebands;
	
	std::string kfacHist;

	int nMasses;

	bool doTraining;
	bool splitSignalSample;
	//int nMassPt;
	std::string names[9];
        std::string BDTnames[9];
	double masses[9];

	std::string mvaWeightsFolder;

protected:
	std::vector<BaseSmearer *> photonSmearers_;
	std::vector<BaseSmearer *> systPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> diPhotonSmearers_;
	std::vector<BaseDiPhotonSmearer *> systDiPhotonSmearers_;
	std::vector<BaseGenLevelSmearer *> genLevelSmearers_;
	std::vector<BaseGenLevelSmearer *> systGenLevelSmearers_;
	
	EnergySmearer /* *eScaleSmearer,*/ *eResolSmearer ; // moved to PhotonAnalysis GF 
	EfficiencySmearer *idEffSmearer, *r9Smearer;
	DiPhoEfficiencySmearer *vtxEffSmearer, *triggerEffSmearer;
	KFactorSmearer * kFactorSmearer;
	
	std::string name_;
	float nevents, sumwei, sumaccept, sumsmear, sumev; 
	
	int nCategories_;
	int nPhotonCategories_;
	int diPhoCounter_;
	// Vertex analysis
	HggVertexAnalyzer vtxAna_;
	HggVertexFromConversions vtxConv_;
	
	// RooStuff
	RooContainer *rooContainer;

	ofstream eventListText;

	//MVA variables
	float _log_H_pt;     
	float _H_ptOverM;     
	float _H_eta;        
	float _d_phi;        
	float _max_eta;      
	float _min_r9;       
	float _pho1_eta;    
	float _pho2_eta;     
	float _pho1_ptOverM; 
	float _pho2_ptOverM; 
	float _sigmaMOverM;
	float _deltaMOverM; 
	float _deltaMOverSigmaM; 
	float _mgg;          
	float _pho1_phi;     
	float _pho1_pt;      
	float _pho1_r9;      
	float _pho2_phi;     
	float _pho2_pt;      
	float _pho2_r9;      
	float _H_pt;         
	float _Ht;           
	float _d_eta;        
	float _mod_d_eta;    
	float _cos_theta_star;
	float _wt;           
	int _cat;           
	int _sideband;           

	//vector<double> weights;
	TFile *kfacFile;
	
    	TMVA::Reader * tmvaReader_;
	TTree * signalTrainTree_[12];
	TTree * signalTestTree_[12];

	TTree * backgroundTrainTree_[12];
	TTree * backgroundTestTree_[12];

	TFile * mvaFile_;
	
};

#endif
