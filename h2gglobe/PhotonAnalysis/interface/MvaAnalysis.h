#ifndef __MVAANALYSIS__
#define __MVAANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "PhotonFix.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
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
	virtual void Analysis(LoopAll&, Int_t);

	double GetPhotonResolution(LoopAll &,int);
	double GetAngleResolutionCorrVtx(LoopAll &,int, int, int);
	double GetAngleResolutionWrongVtx(LoopAll &, int, int, int);
	double PropagateDz(LoopAll &, int, int, int, double);

	double SecH(double);
	double TanH(double);
	
	// Options
	bool reRunCiC;
	float leadEtCut;
	float subleadEtCut;
	std::string efficiencyFile;
	
	// EnergySmearer::energySmearingParameters eSmearPars; // gone to PhotonAnalysis GF
	EfficiencySmearer::efficiencySmearingParameters effSmearPars;
	DiPhoEfficiencySmearer::diPhoEfficiencySmearingParameters diPhoEffSmearPars;

	double GetDifferentialKfactor(double, int);
	bool  doMCSmearing;
	bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
	bool  doEscaleSmear, doEresolSmear, doPhotonIdEffSmear, doVtxEffSmear, doR9Smear, doTriggerEffSmear, doKFactorSmear;
	float systRange;
	int   nSystSteps;   
	int   nEtaCategories, nR9Categories, nPtCategories;
	float massMin, massMax;
	int nDataBins;	
	
	std::string kfacHist;

	TH1D *thm110,*thm120,*thm130,*thm140;
	int nMasses;

	bool doTraining;
	//int nMassPt;
	std::string names[6];
	double masses[6];

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
	float _delta_MOverM; 
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

	//vector<double> weights;
	TFile *kfacFile;
	
    	TMVA::Reader * tmvaReader_;
//        TMVA::Factory* tmvaFactory_[5] ;
	TTree * signalTree_[6];
	TTree * backgroundTree_[6];
	TFile * mvaFile_;
	
};

#endif
