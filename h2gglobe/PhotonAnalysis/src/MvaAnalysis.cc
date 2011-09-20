#include "../interface/MvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::MvaAnalysis()  : 
    name_("MvaAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
    nMasses  = 12;
    signalRegionWidth = 0.07;
}

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::~MvaAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Term(LoopAll& l) 
{

	if (doTraining){
        for (int i = 0; i<nMasses;i++){
            //cout<<"~~Mass = "<<names[i]<<endl;
            //cout<<"test cd"<<endl;
            mvaFile_->cd();
            if (0<signalTree_[i]->GetEntries()    ){
                 //cout<<"test signal: "<<signalTree_[i]->GetEntries()<<endl;
                 signalTree_[i]->Write(("sig"+names[i]).c_str());
            }
            if (0<backgroundTree_[i]->GetEntries()){
                 //cout<<"test background: "<<backgroundTree_[i]->GetEntries()<<endl;
                 backgroundTree_[i]->Write(("bkg"+names[i]).c_str());
            }
        }
        //cout<<"test close"<<endl;
        mvaFile_->Close();
    }
    else{
        float mass_low = masses[2]*(1-signalRegionWidth)/(1+signalRegionWidth);
        //NOTE THE UGLY FIX HERE TO IGNORE THE POINTS 105 and 110
        float mass_high = masses[9]*(1+signalRegionWidth)/(1-signalRegionWidth);
        float mass_boundaries[2];
        mass_boundaries[0] = mass_low*(1-signalRegionWidth);
        mass_boundaries[1] = mass_high*(1+signalRegionWidth);

        for (int i = 2; i<nMasses;i++){
            if (i==8) continue;
            // define hypothesis masses for the sidebands
            float mass_hypothesis = masses[i];
            float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
            float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
            // define the sidebands
            float sideband_boundaries[4];
            sideband_boundaries[0] = 95 > mass_hypothesis_low*(1-signalRegionWidth) ?  95 : mass_hypothesis_low*(1-signalRegionWidth); // only go as low as 95!
            sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
            sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
            sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth); 

            l.rooContainer->FitToData("data_pow_model"+names[i], "data_mass"+names[i],95,sideband_boundaries[1],sideband_boundaries[2],185);
            std::vector<std::pair<double,double> > N_sig = l.rooContainer->GetFitNormalisationsAndErrors("data_pow_model"+names[i],
                                        "data_mass"+names[i],sideband_boundaries[1],sideband_boundaries[2],true);
            
	    cout << mass_hypothesis <<" - pow Normalisation - " << N_sig[0].first << "+/-" << N_sig[0].second << std::endl;
            l.rooContainer->AddNormalisationSystematics("bkg_norm"+names[i],N_sig, 1);

            // Calculate weights to apply to the sidebands
            std::vector<double> wt_low;
            std::vector<double> wt_high;
            for (int i_cat = 0; i_cat<N_sig.size();i_cat++){
                cout<<"N_sig = "<<N_sig[i_cat].first<<endl;
                wt_low.push_back( 0.5*N_sig[i_cat].first);///N_low[i_cat]);    
                wt_high.push_back(0.5*N_sig[i_cat].first);///N_high[i_cat]);    
            }
            // if scale = true, the wt is a scale applied ot the histograms other
            // wise it is an absolute normalisation to be applied
            bool scale = false;//true;
            l.rooContainer->SumBinnedDatasets("data_BDT_sideband_ada"+names[i], "data_low_BDT_ada"+names[i],
                                              "data_high_BDT_ada"+names[i], wt_low, wt_high, scale);
            l.rooContainer->SumBinnedDatasets("data_BDT_sideband_grad"+names[i], "data_low_BDT_grad"+names[i],
                                              "data_high_BDT_grad"+names[i], wt_low, wt_high, scale);

	    std::vector <std::vector<double> > optimizedGradBins; // Need to ompimize binning on MC!
	    optimizedGradBins =  l.rooContainer->OptimizedBinning("bkg_BDT_grad"+names[i],50,false,true,-1);

	    l.rooContainer->RebinBinnedDataset("bkg_grad"+names[i],"data_BDT_sideband_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_high_grad"+names[i],"data_high_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_low_grad"+names[i],"data_low_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_grad"+names[i],"bkg_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_high_grad"+names[i],"bkg_high_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_low_grad"+names[i],"bkg_low_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("data_grad"+names[i],"data_BDT_grad"+names[i],optimizedGradBins,false);
	    l.rooContainer->RebinBinnedDataset("sig_grad"+names[i],"sig_BDT_grad"+names[i],optimizedGradBins,true);

	    std::vector<std::vector <double> > optimizedAdaBins =  l.rooContainer->OptimizedBinning("bkg_BDT_ada"+names[i],50,false,true,-1);
	    l.rooContainer->RebinBinnedDataset("bkg_ada"+names[i],"data_BDT_sideband_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_high_ada"+names[i],"data_high_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_low_ada"+names[i],"data_low_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_ada"+names[i],"bkg_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_high_ada"+names[i],"bkg_high_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("bkg_mc_low_ada"+names[i],"bkg_low_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("data_ada"+names[i],"data_BDT_ada"+names[i],optimizedAdaBins,false);
	    l.rooContainer->RebinBinnedDataset("sig_ada"+names[i],"sig_BDT_ada"+names[i],optimizedAdaBins,true);

            l.rooContainer->WriteDataCard((std::string) l.histFileName+"_ada"+names[i],"data_BDT_ada"+names[i],"sig_BDT_ada"+names[i],"data_BDT_sideband_ada"+names[i]);
            l.rooContainer->WriteDataCard((std::string) l.histFileName+"_grad"+names[i],"data_BDT_grad"+names[i],"sig_BDT_grad"+names[i],"data_BDT_sideband_grad"+names[i]);
        }
    }
	PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
	cout << "InitRealMvaAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;
    
    names[0]="_105";
    masses[0] = 105.;

    names[1]="_110";
    masses[1] = 110.;

    names[2]="_115";
    masses[2] = 115.;

    names[3]="_120";
    masses[3] = 120.;

    names[4]="_125";
    masses[4] = 125.;

    names[5]="_130";
    masses[5] = 130.;

    names[6]="_135";
    masses[6] = 135.;

    names[7]="_140";
    masses[7] = 140.;

    names[8]="_145";
    masses[8] = 145.;

    names[9]="_150";
    masses[9] = 150.;

    names[10]="_121";
    masses[10] = 121.;

    names[11]="_123";
    masses[11] = 123.;

    std::string outputfilename = (std::string) l.histFileName;
    eventListText.open(Form("%s_ascii_events.txt",outputfilename.c_str()));
    //
    // These parameters are set in the configuration file
    std::cout
	<< "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "MvaAnalysis " << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< "leadEtCut "<< leadEtCut << "\n"
	<< "subleadEtCut "<< subleadEtCut << "\n"
	<< "doTriggerSelection "<< doTriggerSelection << "\n"
	<< "nEtaCategories "<< nEtaCategories << "\n"
	<< "nR9Categories "<< nR9Categories << "\n"		
	<< "nPtCategories "<< nPtCategories << "\n"		
	<< "doEscaleSyst "<< doEscaleSyst << "\n"
	<< "doEresolSyst "<< doEresolSyst << "\n"
	<< "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
	<< "efficiencyFile " << efficiencyFile << "\n"
	<< "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
	<< "doR9Syst "<< doR9Syst << "\n"
	<< "doVtxEffSyst "<< doVtxEffSyst << "\n"
	<< "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
	<< "doKFactorSyst "<< doKFactorSyst << "\n"
	<< "-------------------------------------------------------------------------------------- \n"
	<< std::endl;

    // avoid recalculated the CIC ID every time
    l.runCiC = reRunCiC;
    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
	l.histoContainer[ind].setScale(1.);
    }
    
    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;
    

    effSmearPars.categoryType = "2CatR9_EBEE";
    effSmearPars.n_categories = 4;
    effSmearPars.efficiency_file = efficiencyFile;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFile;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
        photonSmearers_.push_back(eScaleSmearer);
    }
    if( doEresolSmear ) {
	// energy resolution smearing
	std::cerr << __LINE__ << std::endl; 
	eResolSmearer = new EnergySmearer( eSmearPars );
	eResolSmearer->name("E_res");
	eResolSmearer->doEnergy(false);
	eResolSmearer->scaleOrSmear(false);
	photonSmearers_.push_back(eResolSmearer);
    }
    if( doPhotonIdEffSmear ) {
	// photon ID efficiency 
	std::cerr << __LINE__ << std::endl; 
	idEffSmearer = new EfficiencySmearer( effSmearPars );
	idEffSmearer->name("idEff");
	idEffSmearer->setEffName("ratioTP");
	idEffSmearer->init();
	idEffSmearer->doPhoId(true);
	photonSmearers_.push_back(idEffSmearer);
    }
    if( doR9Smear ) {
	// R9 re-weighting
	r9Smearer = new EfficiencySmearer( effSmearPars );
	r9Smearer->name("r9Eff");
	r9Smearer->setEffName("ratioR9");
	r9Smearer->init();
	r9Smearer->doR9(true);
	photonSmearers_.push_back(r9Smearer);
    }
    if( doVtxEffSmear ) {
	// Vertex ID
	std::cerr << __LINE__ << std::endl; 
	vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
	vtxEffSmearer->name("vtxEff");
	vtxEffSmearer->setEffName("ratioVertex");
	vtxEffSmearer->doVtxEff(true);
	vtxEffSmearer->init();
	diPhotonSmearers_.push_back(vtxEffSmearer);
    }
    if( doTriggerEffSmear ) {
	// trigger efficiency
	std::cerr << __LINE__ << std::endl; 
	triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
	triggerEffSmearer->name("triggerEff");
	triggerEffSmearer->setEffName("effL1HLT");
	triggerEffSmearer->doVtxEff(false);
	triggerEffSmearer->init();
	diPhotonSmearers_.push_back(triggerEffSmearer);
    }
    if(doKFactorSmear) {
	// kFactor efficiency
	std::cerr << __LINE__ << std::endl; 
	kFactorSmearer = new KFactorSmearer( kfacHist );
	kFactorSmearer->name("kFactor");
	kFactorSmearer->init();
	genLevelSmearers_.push_back(kFactorSmearer);
    }


    // RooContainer stuff
    // FIXME move these params to config file
    //l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->SetNCategories(1);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    // RooContainer does not support steps different from 1 sigma
    //assert( ((float)nSystSteps) == systRange );
    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
	std::vector<std::string> sys(1,eCorrSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
	systPhotonSmearers_.push_back( eScaleSmearer );
	std::vector<std::string> sys(1,eScaleSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
	systPhotonSmearers_.push_back( eResolSmearer );
	std::vector<std::string> sys(1,eResolSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
	systPhotonSmearers_.push_back( idEffSmearer );
	std::vector<std::string> sys(1,idEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
	systPhotonSmearers_.push_back( r9Smearer );
	std::vector<std::string> sys(1,r9Smearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
	systDiPhotonSmearers_.push_back( vtxEffSmearer );
	std::vector<std::string> sys(1,vtxEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
	systDiPhotonSmearers_.push_back( triggerEffSmearer );
	std::vector<std::string> sys(1,triggerEffSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
	systGenLevelSmearers_.push_back(kFactorSmearer);
	std::vector<std::string> sys(1,kFactorSmearer->name());
	std::vector<int> sys_t(1,-1);	// -1 for signal, 1 for background 0 for both
	l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
	
    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.06,1.00);
    // ----------------------------------------------------

    // Create observables for shape-analysis with ranges
    // l.rooContainer->AddObservable("mass" ,100.,150.);
    //mass low is centre of lowest sideband for lowest test mass
    //mass hgih is centre of hgihest sideband for hgihest test mass
    float mass_low = masses[2]*(1-signalRegionWidth)/(1+signalRegionWidth);
    //NOTE THE UGLY FIX HERE TO IGNORE THE POINTS 105 and 110
    float mass_high = masses[9]*(1+signalRegionWidth)/(1-signalRegionWidth);
    float mass_boundaries[2];
    mass_boundaries[0] = mass_low*(1-signalRegionWidth);
    mass_boundaries[1] = mass_high*(1+signalRegionWidth);

    l.rooContainer->AddObservable("CMS_hgg_mass",95,185);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    for (int i = 2; i<nMasses;i++){

    l.rooContainer->AddRealVar("mu"+names[i],-0.01,-0.15,0);
    l.rooContainer->AddRealVar("mu1"+names[i],-0.01,-0.15,0);
    l.rooContainer->AddRealVar("mu2"+names[i],-0.05,-0.15,0);
    l.rooContainer->AddRealVar("mu1a"+names[i],-0.01,-0.15,0);
    l.rooContainer->AddRealVar("mu2a"+names[i],-0.05,-0.15,0);
    l.rooContainer->AddRealVar("mu3a"+names[i],-0.08,-0.15,0);
    l.rooContainer->AddRealVar("f2"+names[i],0.3,0.,1.);
    l.rooContainer->AddRealVar("f2a"+names[i],0.3,0.,1.);
    l.rooContainer->AddRealVar("f3a"+names[i],0.3,0.,1.);
    l.rooContainer->AddRealVar("r1"+names[i],-1,-5.,0);
    l.rooContainer->AddRealVar("f1"+names[i],1.,0.,2.);
    l.rooContainer->AddRealVar("r2"+names[i],-1,-5.,0);
    l.rooContainer->AddRealVar("f2"+names[i],1.,0.,2.);
    l.rooContainer->AddRealVar("pol0"+names[i],1,-1.5,1.5);
    l.rooContainer->AddRealVar("pol1"+names[i],1,-1.5,1.5);
    l.rooContainer->AddRealVar("pol2"+names[i],1,-1.5,1.5);
    //l.rooContainer->AddRealVar("pol1",-0.01,-1.5,1.5);
    l.rooContainer->AddFormulaVar("modpol0"+names[i],"@0*@0","pol0"+names[i]);
    l.rooContainer->AddFormulaVar("modpol1"+names[i],"@0*@0","pol1"+names[i]);

    // Exponential Model
    std::vector<std::string> data_exp_pars(1,"p");	 
    data_exp_pars[0] = "mu"+names[i];
    l.rooContainer->AddGenericPdf("data_pol_model"+names[i], "0","CMS_hgg_mass",data_exp_pars,1);

    // 2nd order polynomial model
    std::vector<std::string> data_pol_pars(2,"p");	 
    data_pol_pars[0] = "modpol0"+names[i];
    data_pol_pars[1] = "modpol1"+names[i];
    l.rooContainer->AddGenericPdf("data_poly_model"+names[i], "0","CMS_hgg_mass",data_pol_pars,72);	

    // Double Exponential model  -- The RooAddPdf Way
    std::vector<std::string> data_exp1_pars(1,"p");	 
    data_exp1_pars[0] = "mu1"+names[i];
    l.rooContainer->AddGenericPdf("data_expa_model"+names[i], "0","CMS_hgg_mass",data_exp1_pars,1,0.49,0.0,0.9);

    std::vector<std::string> data_exp2_pars(1,"p");	 
    data_exp2_pars[0] = "mu2"+names[i];
    l.rooContainer->AddGenericPdf("data_expb_model"+names[i], "0","CMS_hgg_mass",data_exp2_pars,1,0.49,0.0,0.9);

    std::vector<std::string> components_data(2,"c");
    components_data[1]="data_expa_model"+names[i];	
    components_data[0]="data_expb_model"+names[i];	
    l.rooContainer->ComposePdf("data_exp2_model"+names[i],"exp1+exp2",components_data,false);

    // Double Exponential model  -- The Formula Way
    /*
    std::vector<std::string> data_exp1_pars(4,"p");	 
    data_exp1_pars[0] = "mu1";
    data_exp1_pars[1] = "f1";
    data_exp1_pars[2] = "mu2";
    data_exp1_pars[3] = "f2";
    l.rooContainer->AddGenericPdf("data_exp2_model"+names[i], "@2*TMath::Exp(@0*@1) + @4*TMath::Exp(@0*@3)","CMS_hgg_mass",data_exp1_pars,0);
    */


    // Triple Exponential model
    std::vector<std::string> data_expa1_pars(1,"p");	 
    data_expa1_pars[0] = "mu1a"+names[i];
    l.rooContainer->AddGenericPdf("data_expaa_model"+names[i], "0","CMS_hgg_mass",data_expa1_pars,1,0.33,0.0,0.9);

    std::vector<std::string> data_expa2_pars(1,"p");	 
    data_expa2_pars[0] = "mu2a"+names[i];
    l.rooContainer->AddGenericPdf("data_expbb_model"+names[i], "0","CMS_hgg_mass",data_expa2_pars,1,0.29,0.0,0.9);

    std::vector<std::string> data_expa3_pars(1,"p");	 
    data_expa3_pars[0] = "mu3a"+names[i];
    l.rooContainer->AddGenericPdf("data_expcc_model"+names[i], "0","CMS_hgg_mass",data_expa3_pars,1,0.29,0.0,0.9);

    std::vector<std::string> components_3data(3,"c");
    components_3data[2]="data_expaa_model"+names[i];	
    components_3data[1]="data_expbb_model"+names[i];	
    components_3data[0]="data_expcc_model"+names[i];	
    l.rooContainer->ComposePdf("data_exp3_model"+names[i],"exp1+exp2+exp3",components_3data,false);

    // Laurent series
    std::vector<std::string> data_lau_pars(2,"p");	 
    data_lau_pars[0] = "f1"+names[i];
    data_lau_pars[1] = "f2"+names[i];
    l.rooContainer->AddGenericPdf("data_lau_model"+names[i], "(@1/(@0*@0*@0*@0)) + (@2/(@0*@0*@0*@0*@0)) ","CMS_hgg_mass",data_lau_pars,0);

    // Power law
    std::vector<std::string> data_pow_pars(1,"p");	 
    data_pow_pars[0] = "r1"+names[i];
    l.rooContainer->AddGenericPdf("data_pow_model"+names[i], "TMath::Power(@0,@1)","CMS_hgg_mass",data_pow_pars,0);

    }
        
	
	if (doTraining){
	    TString outfileName( "TMVA_input_" + (std::string) l.histFileName);
        mvaFile_ = TFile::Open( outfileName, "RECREATE" );
        mvaFile_->cd();
        for (int i = 0; i<nMasses;i++){
            signalTree_[i]          = new TTree(("sig"+names[i]).c_str(),
                                            ("SignalTree"+names[i]).c_str());
            TBranch *s_log_H_pt     = signalTree_[i]->Branch("log_H_pt", &_log_H_pt , "log_H_pt/F");;
            TBranch *s_H_ptOverM    = signalTree_[i]->Branch("H_ptOverM", &_H_ptOverM , "_H_ptOverM/F");;
            TBranch *s_H_eta        = signalTree_[i]->Branch("H_eta", &_H_eta , "H_eta/F");;
            TBranch *s_d_phi        = signalTree_[i]->Branch("d_phi", &_d_phi,"d_phi/F");
            TBranch *s_max_eta      = signalTree_[i]->Branch("max_eta", &_max_eta , "max_eta/F");;
            TBranch *s_min_r9       = signalTree_[i]->Branch("min_r9", &_min_r9 , "min_r9/F");;
            TBranch *s_pho1_eta     = signalTree_[i]->Branch("pho1_eta", &_pho1_eta , "pho1_eta/F");
            TBranch *s_pho2_eta     = signalTree_[i]->Branch("pho2_eta", &_pho2_eta , "pho2_eta/F");
            TBranch *s_pho1_ptOverM = signalTree_[i]->Branch("pho1_ptOverM", &_pho1_ptOverM , "pho1_ptOverM/F");;
            TBranch *s_pho2_ptOverM = signalTree_[i]->Branch("pho2_ptOverM", &_pho2_ptOverM , "pho2_ptOverM/F");;
            TBranch *s_deltaMOverM = signalTree_[i]->Branch("deltaMOverM", &_deltaMOverM,"deltaMOverM/F");
            TBranch *s_deltaMOverSigmaM = signalTree_[i]->Branch("deltaMOverSigmaM", &_deltaMOverSigmaM,"deltaMOverSigmaM/F");
            TBranch *s_sigmaMOverM = signalTree_[i]->Branch("sigmaMOverM", &_sigmaMOverM,"sigmaMOverM/F");
            TBranch *s_mgg          = signalTree_[i]->Branch("mgg", &_mgg, "mgg/F");
            TBranch *s_pho1_phi     = signalTree_[i]->Branch("pho1_phi", &_pho1_phi , "pho1_phi/F");
            TBranch *s_pho1_pt      = signalTree_[i]->Branch("pho1_pt", &_pho1_pt , "pho1_pt/F");;
            TBranch *s_pho1_r9      = signalTree_[i]->Branch("pho1_r9", &_pho1_r9 , "_pho1_r9/F");;
            TBranch *s_pho2_phi     = signalTree_[i]->Branch("pho2_phi", &_pho2_phi , "pho2_phi/F");
            TBranch *s_pho2_pt      = signalTree_[i]->Branch("pho2_pt", &_pho2_pt , "pho2_pt/F");;
            TBranch *s_pho2_r9      = signalTree_[i]->Branch("pho2_r9", &_pho2_r9 , "_pho2_r9/F");;
            TBranch *s_H_pt         = signalTree_[i]->Branch("H_pt", &_H_pt , "H_pt/F");;
            TBranch *s_Ht           = signalTree_[i]->Branch("Ht", &_Ht , "Ht/F");;
            TBranch *s_d_eta        = signalTree_[i]->Branch("d_eta", &_d_eta,"d_eta/F");
            TBranch *s_mod_d_eta    = signalTree_[i]->Branch("mod_d_eta", &_mod_d_eta,"mod_d_eta/F");
            TBranch *s_cos_theta_star= signalTree_[i]->Branch("cos_theta_star", &_cos_theta_star , "cos_theta_star/F");;
            TBranch *s_wt           = signalTree_[i]->Branch("wt", &_wt, "wt/F");

            backgroundTree_[i]      = new TTree(("bkg"+names[i]).c_str(),
                                                ("BackgroundTree"+names[i]).c_str());
            TBranch *b_log_H_pt     = backgroundTree_[i]->Branch("log_H_pt", &_log_H_pt , "log_H_pt/F");;
            TBranch *b_H_ptOverM    = backgroundTree_[i]->Branch("H_ptOverM", &_H_ptOverM , "_H_ptOverM/F");;
            TBranch *b_H_eta        = backgroundTree_[i]->Branch("H_eta", &_H_eta , "H_eta/F");;
            TBranch *b_d_phi        = backgroundTree_[i]->Branch("d_phi", &_d_phi,"d_phi/F");
            TBranch *b_max_eta      = backgroundTree_[i]->Branch("max_eta", &_max_eta , "max_eta/F");;
            TBranch *b_min_r9       = backgroundTree_[i]->Branch("min_r9", &_min_r9 , "min_r9/F");;
            TBranch *b_pho1_eta     = backgroundTree_[i]->Branch("pho1_eta", &_pho1_eta , "pho1_eta/F");
            TBranch *b_pho2_eta     = backgroundTree_[i]->Branch("pho2_eta", &_pho2_eta , "pho2_eta/F");
            TBranch *b_pho1_ptOverM = backgroundTree_[i]->Branch("pho1_ptOverM", &_pho1_ptOverM , "pho1_ptOverM/F");;
            TBranch *b_pho2_ptOverM = backgroundTree_[i]->Branch("pho2_ptOverM", &_pho2_ptOverM , "pho2_ptOverM/F");;
            TBranch *b_deltaMOverM = backgroundTree_[i]->Branch("deltaMOverM", &_deltaMOverM,"deltaMOverM/F");
            TBranch *b_deltaMOverSigmaM = backgroundTree_[i]->Branch("deltaMOverSigmaM", &_deltaMOverSigmaM,"deltaMOverSigmaM/F");
            TBranch *b_sigmaMOverM = backgroundTree_[i]->Branch("sigmaMOverM", &_sigmaMOverM,"sigmaMOverM/F");
            TBranch *b_mgg          = backgroundTree_[i]->Branch("mgg", &_mgg, "mgg/F");
            TBranch *b_pho1_phi     = backgroundTree_[i]->Branch("pho1_phi", &_pho1_phi , "pho1_phi/F");
            TBranch *b_pho1_pt      = backgroundTree_[i]->Branch("pho1_pt", &_pho1_pt , "pho1_pt/F");;
            TBranch *b_pho1_r9      = backgroundTree_[i]->Branch("pho1_r9", &_pho1_r9 , "_pho1_r9/F");;
            TBranch *b_pho2_phi     = backgroundTree_[i]->Branch("pho2_phi", &_pho2_phi , "pho2_phi/F");
            TBranch *b_pho2_pt      = backgroundTree_[i]->Branch("pho2_pt", &_pho2_pt , "pho2_pt/F");;
            TBranch *b_pho2_r9      = backgroundTree_[i]->Branch("pho2_r9", &_pho2_r9 , "_pho2_r9/F");;
            TBranch *b_H_pt         = backgroundTree_[i]->Branch("H_pt", &_H_pt , "H_pt/F");;
            TBranch *b_Ht           = backgroundTree_[i]->Branch("Ht", &_Ht , "Ht/F");;
            TBranch *b_d_eta        = backgroundTree_[i]->Branch("d_eta", &_d_eta,"d_eta/F");
            TBranch *b_mod_d_eta    = backgroundTree_[i]->Branch("mod_d_eta", &_mod_d_eta,"mod_d_eta/F");
            TBranch *b_cos_theta_star= backgroundTree_[i]->Branch("cos_theta_star", &_cos_theta_star , "cos_theta_star/F");;
            TBranch *b_wt           = backgroundTree_[i]->Branch("wt", &_wt, "wt/F");
        }
    }
    else{

        l.rooContainer->AddObservable("BDT" ,-1.,1.);

        //Set up TMVA reader
        tmvaReader_= new TMVA::Reader();
 		tmvaReader_->AddVariable("H_ptOverM", &_H_ptOverM);
 		tmvaReader_->AddVariable("H_eta", &_H_eta);
 		tmvaReader_->AddVariable("d_phi", &_d_phi);
 		tmvaReader_->AddVariable("max_eta", &_max_eta);
 		tmvaReader_->AddVariable("min_r9", &_min_r9);
 		tmvaReader_->AddVariable("pho1_eta", &_pho1_eta);
 		tmvaReader_->AddVariable("pho2_eta", &_pho2_eta);
 		tmvaReader_->AddVariable("pho1_ptOverM", &_pho1_ptOverM);
 		tmvaReader_->AddVariable("pho2_ptOverM", &_pho2_ptOverM);
 	    tmvaReader_->AddVariable("deltaMOverM", &_deltaMOverM);
 	    //tmvaReader_->AddVariable("deltaMOverSigmaM", &_deltaMOverSigmaM);
		// TEST TO SEE DIFFERENCE WITHOUT SIGMA M
 		tmvaReader_->AddVariable("sigmaMOverM", &_sigmaMOverM);

        for (int i = 2; i<nMasses;i++){  // We are ignoring masses 105 and 110 for now
            if (i==8) continue;//Not available yet
	    int nBDTbins = 5000;
            //Adaptive Boost
            l.rooContainer->CreateDataSet("BDT","data_low_BDT_ada"+names[i]  ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","data_BDT_ada"+names[i]	     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","data_high_BDT_ada"+names[i] ,nBDTbins);

            l.rooContainer->CreateDataSet("BDT","bkg_low_BDT_ada"+names[i]   ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","bkg_BDT_ada"+names[i]       ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","bkg_high_BDT_ada"+names[i]  ,nBDTbins);

            l.rooContainer->CreateDataSet("BDT","sig_BDT_ada"+names[i]       ,nBDTbins);    

            //Gradiant Boost
            l.rooContainer->CreateDataSet("BDT","data_low_BDT_grad"+names[i] ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","data_BDT_grad"+names[i]     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","data_high_BDT_grad"+names[i],nBDTbins);

            l.rooContainer->CreateDataSet("BDT","bkg_low_BDT_grad"+names[i]  ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","bkg_BDT_grad"+names[i]      ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT","bkg_high_BDT_grad"+names[i] ,nBDTbins);

            l.rooContainer->CreateDataSet("BDT","sig_BDT_grad"+names[i]      ,nBDTbins);    

            //Invariant Mass Spectra
            l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass"+names[i],nDataBins); // (100,110,150) -> for a window, else full obs range is taken 
            l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass"+names[i] ,nDataBins);    	 
		
            // Make the signal Systematic Sets
            l.rooContainer->MakeSystematics("BDT","sig_BDT_grad"+names[i] ,-1)	;
            l.rooContainer->MakeSystematics("BDT","sig_BDT_ada"+names[i]  ,-1)	;

            //TMVA Reader
            tmvaReader_->BookMVA("BDT_ada" +names[i],mvaWeightsFolder+"/TMVAClassification_BDT_ada" +names[i]+".weights.xml");
            tmvaReader_->BookMVA("BDT_grad"+names[i],mvaWeightsFolder+"/TMVAClassification_BDT_grad"+names[i]+".weights.xml");
        }
    }

    //PhotonFix::initialise("Nominal");
    /* -----------------------------------------------------------------------------------------
       KFactors Reweighting
       ------------------------------------------------------------------------------------------- */
    if(PADEBUG) 
	cout << "InitRealMvaAnalysis END"<<endl;
	
    // FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
    if(PADEBUG) 
	cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;
    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    //PU reweighting
    unsigned int n_pu = l.pu_n;
    if ( cur_type !=0 && puHist != "") {
	bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
	if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
	    std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
	}
	std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
	if(n_pu<puweights.size()){
	    weight *= puweights[n_pu];
	    sumwei+=puweights[n_pu]; 
	}    
	else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	    cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
	}
    }
    
    assert( weight >= 0. );  
    l.FillCounter( "PUWeighted", weight );
    
    if( jentry % 10000 ==  0 ) {
	    std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
		      << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
		      <<  sumaccept / sumev << " " << sumsmear / sumaccept
		      << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){            // if background sample, gP4 remains 4vect(0)
	for (int gi=0;gi<l.gp_n;gi++){
	    if (l.gp_pdgid[gi]==25){
		gP4 = *((TLorentzVector*)l.gp_p4->At(gi));
		gPT = gP4.Pt();
		break;
	    }
	}
    }

    // ------------------------------------------------------------
    // smear all of the photons!
    std::pair<int,int> diphoton_index;
   
    // do gen-level dependent first (e.g. k-factor); only for signal
    double genLevWeight=1; 
    if(cur_type<0){
	for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
	    float genWeight=1;
	    (*si)->smearEvent( genWeight,gP4, l.pu_n, cur_type, 0. );
	    if( genWeight < 0. ) {
		std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		assert(0);
	    }
	    genLevWeight*=genWeight;
	}
    }

    // Nominal smearing
    std::vector<float> smeared_pho_energy(l.pho_n,0.); 
    std::vector<float> smeared_pho_r9(l.pho_n,0.); 
    std::vector<float> smeared_pho_weight(l.pho_n,1.);
   
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
	std::vector<std::vector<bool> > p;
	PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
				    // *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
				    ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_residCorrEnergy[ipho],
				    l.pho_isEB[ipho], l.pho_r9[ipho],
				    l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_) );
	float pweight = 1.;
	// smear MC. But apply energy shift to data 
	if( cur_type != 0 && doMCSmearing ) { // if it's MC
	    for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
		float sweight = 1.;
		(*si)->smearPhoton(phoInfo,sweight,l.run,0.);	   
		if( sweight < 0. ) {
			std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
			assert(0);
		}
		pweight *= sweight;
	    }
	} else if(cur_type == 0 ) {          // if it's data
	    if (doEcorrectionSmear){
	      float sweight = 1.; 
	      eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.);     // This Smearer is the same as for MC so can just re-use it
	      pweight *= sweight;	
	    }
 	    if (doEscaleSmear){
	      float sweight = 1.; 
	      eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
	      pweight *= sweight;
	    }
	}

	smeared_pho_energy[ipho] = phoInfo.energy();
	smeared_pho_r9[ipho] = phoInfo.r9();
	smeared_pho_weight[ipho] = pweight;
    }
   
    sumev += weight;
    // FIXME pass smeared R9
    int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
    /// std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id > -1 ) {

	diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
    	// bring all the weights together: lumi & Xsection, single gammas, pt kfactor
	float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

	l.countersred[diPhoCounter_]++;

	TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

	bool CorrectVertex;
	// FIXME pass smeared R9
	int category = 0;//l.gge Society says thDiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,nPtCategories);
	int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
	if( cur_type != 0 && doMCSmearing ) {
	    float pth = Higgs.Pt();
	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
		float rewei=1.;
		(*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
		if( rewei < 0. ) {
		    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
		    assert(0);
		}
		evweight *= rewei;
	    }
	    CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
	}
	float mass    = Higgs.M();
	float ptHiggs = Higgs.Pt();
      
	assert( evweight >= 0. ); 

/*
	double lead_photonResolution = GetPhotonResolution(l,diphoton_index.first);
	double sublead_photonResolution = GetPhotonResolution(l,diphoton_index.second);
	double angle_resolution = GetAngleResolutionCorrVtx(l,diphoton_index.first,diphoton_index.second,diphoton_id);
*/
	
	
	// returns massResolution assuming no error on vertex (i.e. just from Paul's stuff)
	//double massResolutionEonly = 0.5*mass*TMath::Sqrt((lead_photonResolution*lead_photonResolution)/(lead_p4.E()*lead_p4.E())
	//					 +(sublead_photonResolution*sublead_photonResolution)/(sublead_p4.E()*sublead_p4.E()));
  
	//double alpha = lead_p4.Angle(sublead_p4.Vect());

	// returns massResolution with vertex error (at the moment only assumes correct vertex)
	//double massResolution = 0.5*mass*TMath::Sqrt((lead_photonResolution*lead_photonResolution)/(lead_p4.E()*lead_p4.E())
	//	+(sublead_photonResolution*sublead_photonResolution)/(sublead_p4.E()*sublead_p4.E())
	//	+((angle_resolution*angle_resolution)*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))*(TMath::Sin(alpha)/(1.-TMath::Cos(alpha)))));

	// Mass Resolution of the Event
	massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

	double massResolution = massResolutionCalculator->massResolution();

        if (doTraining){
            if (cur_type > 0 ){// Background 
                for (int i = 0; i<nMasses;i++) {
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,masses[i],evweight);
                    backgroundTree_[i]->Fill();
                } 
            }
            else { //Signal 
                int i0 = SignalType(cur_type); 
                if (i0<0){
                    cout<<"CAN'T TRAIN ON DATA!\n";
                    return;
                }
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,masses[i0],evweight);
                signalTree_[i0]->Fill();
            }
        }
        else{
            // Iterate over each mass point. 
            for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                if (i==8) continue;
			l.rooContainer->InputDataPoint("data_mass"+names[i],category,mass,evweight);
		}

                // define hypothesis masses for the sidebands
                float mass_hypothesis = masses[i];
                float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                // define the sidebands
                float sideband_boundaries[4];
                sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);

                //Signal Window
                if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2]){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);

                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );

                    if (cur_type == 0 ){//data

			
                        l.FillHist("deltaMOverM"+names[i],0, _deltaMOverM, evweight);
                        l.FillHist("deltaMOverSigmaM"+names[i],0, _deltaMOverSigmaM, evweight);
                        l.FillHist("sigmaMOverM"+names[i],0, _sigmaMOverM, evweight);
                        l.FillHist("pho1_eta"+names[i],0, _pho1_eta, evweight);
                        l.FillHist("pho2_eta"+names[i],0, _pho2_eta, evweight);
                        l.FillHist("pho1_ptOverM"+names[i],0, _pho1_ptOverM, evweight);
                        l.FillHist("pho2_ptOverM"+names[i],0, _pho2_ptOverM, evweight);
                        l.FillHist("log_H_pt"+names[i],0, _log_H_pt, evweight);
                        l.FillHist("H_eta"+names[i],0, _H_eta, evweight);
                        l.FillHist("d_phi"+names[i],0, _d_phi, evweight);
                        l.FillHist("max_eta"+names[i],0, _max_eta, evweight);
                        l.FillHist("min_r9"+names[i],0, _min_r9, evweight);

                        l.rooContainer->InputDataPoint("data_BDT_ada"+names[i],category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("data_BDT_grad"+names[i] ,category,bdt_grad,evweight);

                        l.FillHist("BDT_ada"+names[i],0, bdt_ada, evweight);
                        l.FillHist("BDT_grad"+names[i],0, bdt_grad, evweight);
                    }
                    if (cur_type > 0 ){// background MC
                        l.rooContainer->InputDataPoint("bkg_BDT_ada"+names[i],category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("bkg_BDT_grad"+names[i] ,category,bdt_grad,evweight);
                    }
                    if (cur_type < 0){// signal MC
                        // Fill if the current type of MC matches the
                        // current iteration over the masses 
                        if (SignalType(cur_type)==i){
                            l.rooContainer->InputDataPoint("sig_BDT_ada"+names[i],category,bdt_ada,evweight);
                            l.rooContainer->InputDataPoint("sig_BDT_grad"+names[i] ,category,bdt_grad,evweight);
                        }
                    }
                }
                //Lower Window
                else if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_low,evweight);
                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );
                    if (cur_type == 0 ){//data
                         l.rooContainer->InputDataPoint("data_low_BDT_ada"+names[i],category,bdt_ada,evweight);
                         l.rooContainer->InputDataPoint("data_low_BDT_grad"+names[i] ,category,bdt_grad,evweight);
                    }
                    if (cur_type > 0 ){// background MC
                        l.rooContainer->InputDataPoint("bkg_low_BDT_ada"+names[i] ,category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("bkg_low_BDT_grad"+names[i],category,bdt_grad,evweight);
                    }
                }
                //Upper Window
                else if( mass>sideband_boundaries[2] && mass<sideband_boundaries[3]){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis_high,evweight);
                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );
                        l.rooContainer->InputDataPoint("data_high_BDT_ada"+names[i],category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("data_high_BDT_grad"+names[i] ,category,bdt_grad,evweight);
                    }
                    if (cur_type > 0 ){// background MC
                        l.rooContainer->InputDataPoint("bkg_high_BDT_ada"+names[i] ,category,bdt_ada,evweight);
                        l.rooContainer->InputDataPoint("bkg_high_BDT_grad"+names[i],category,bdt_grad,evweight);
                    }
                }
            }
        }
    //}
	l.FillCounter( "Accepted", weight );
	l.FillCounter( "Smeared", evweight );
	sumaccept += weight;
 	sumsmear += evweight;
    }
    if(PADEBUG) 
	cout<<"myFillHistRed END"<<endl;


    // Now do some MC Systematic studies - For MVA Analysis, assume that we only care about Signal Systematics Here, ie make the check cur_type < 0 to increase speed
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
    if( cur_type < 0 && doMCSmearing && !doTraining){
	// fill steps for syst uncertainty study
	float systStep = systRange / (float)nSystSteps;
	// di-photon smearers systematics
	if (diphoton_id > -1 ) {
	       
	    float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
	    float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
	    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
	    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
	 
	    for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
		std::vector<double> bdt_grad_errors;
		std::vector<double> bdt_ada_errors;
		std::vector<double> weights;
		std::vector<int>    categories;
	   
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
	     
		    int category = 0;  // Category always 0 fro MVA analysis
		    double genLevWeightSyst=1; 
	     
		    for(std::vector<BaseGenLevelSmearer *>::iterator sj=genLevelSmearers_.begin(); sj!= genLevelSmearers_.end(); ++sj ) {
			float swei=1.;
			if( *si == *sj ) { 
			    (*si)->smearEvent(swei, gP4, l.pu_n, cur_type, syst_shift );
			} else {
			    (*sj)->smearEvent(swei, gP4, l.pu_n, cur_type, 0. );
			}
			genLevWeightSyst *= swei;
		    }

		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeightSyst;
		    float mass = Higgs.M();
		    float ptHiggs = Higgs.Pt();

		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();

                   // Iterate over each mass point. 
            for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                if (SignalType(cur_type)!=i) continue;

                // define hypothesis masses for the sidebands
                float mass_hypothesis = masses[i];
                float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                // define the sidebands
                float sideband_boundaries[4];
                sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);

                //Signal Window
                if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2]){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                    float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                    float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );

                    categories.push_back(category);
                    bdt_ada_errors.push_back(bdt_ada);
                    bdt_grad_errors.push_back(bdt_grad);
                    weights.push_back(evweight);

                } else {

                    categories.push_back(-1);
                    bdt_ada_errors.push_back(-100.);
                    bdt_grad_errors.push_back(-100.);
                    weights.push_back(0.);
                }
            }
	    }// end loop on systematics steps

		// Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                // Iterate over each mass point. 
                for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                  if (SignalType(cur_type)==i){
                      l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i],(*si)->name(),categories,bdt_ada_errors,weights);
                      l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i],(*si)->name(),categories,bdt_grad_errors,weights);
                 }
	        }
		// -------------------------------------------------------------------------------------------------------------//
 
	    } // end loop on smearers 
		 

	    for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
		std::vector<double> bdt_grad_errors;
		std::vector<double> bdt_ada_errors;
		std::vector<double> weights;
		std::vector<int> categories;
		       
		for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		    if( syst_shift == 0. ) { continue; } // skip the central value
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
			       
		    // restart with 'fresh' wait for this round of systematics
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;
			       
		    // FIXME pass smeared R9 and di-photon
		    int category = 0; 
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
			float swei=1.;
			float pth = Higgs.Pt();
			if( *si == *sj ) { 
			    (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), syst_shift );
			} else { 
			    (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			}
			evweight *= swei;
		    }
		    float mass = Higgs.M();
		    float ptHiggs = Higgs.Pt();
		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();
                   // Iterate over each mass point. 
                   for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                     if (SignalType(cur_type)!=i) continue;

                     // define hypothesis masses for the sidebands
                     float mass_hypothesis = masses[i];
                     float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                     float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                     // define the sidebands
                     float sideband_boundaries[4];
                     sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                     sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                     sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                     sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);
		
                     //Signal Window
                     if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2]){//Signal mass window cut
                       SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                       float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                       float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );

		       categories.push_back(category);
		       bdt_ada_errors.push_back(bdt_ada);
		       bdt_grad_errors.push_back(bdt_grad);
		       weights.push_back(evweight);

		    } else {

		       categories.push_back(-1);
		       bdt_ada_errors.push_back(-100.);
		       bdt_grad_errors.push_back(-100.);
		       weights.push_back(0.);
		    }
		  }
	        }// end loop on systematics steps

		// Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                // Iterate over each mass point. 
                for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                  if (SignalType(cur_type)==i){
                      l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i],(*si)->name(),categories,bdt_ada_errors,weights);
                      l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i],(*si)->name(),categories,bdt_grad_errors,weights);
                 }
	        }
		// -------------------------------------------------------------------------------------------------------------//
 
	    } // end loop on smearers 
		 
	} // Close If on CiC Selection
	// Now the systematic Steps which can effect the CiC Selection

       
	// loop over the smearers included in the systematics study
	for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
	    std::vector<double> bdt_grad_errors;
	    std::vector<double> bdt_ada_errors;
	    std::vector<double> weights;
	    std::vector<int> categories;
	   
	    // loop over syst shift
	    for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
		if( syst_shift == 0. ) { continue; } // skip the central value
		// smear the photons 
		for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
		    std::vector<std::vector<bool> > p;
		    //std::cout << "GF check: " <<  l.pho_residCorrEnergy[ipho] << "  " << l.pho_residCorrResn[ipho] << std::endl;
		    PhotonReducedInfo phoInfo ( *((TVector3*)l.pho_calopos->At(ipho)), 
						/// *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
						((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), l.pho_residCorrEnergy[ipho],
						l.pho_isEB[ipho], l.pho_r9[ipho],
						l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_));
		  
		    float pweight = 1.;
		    for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
			float sweight = 1.;
			if( *si == *sj ) {
			    // move the smearer under study by syst_shift
			    (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
			} else {
			    // for the other use the nominal points
			    (*sj)->smearPhoton(phoInfo,sweight,l.run,0.);
			}
			pweight *= sweight;
		    }
		    smeared_pho_energy[ipho] = phoInfo.energy();
		    smeared_pho_r9[ipho] = phoInfo.r9();
		    smeared_pho_weight[ipho] = pweight;
		}
	       
		// analyze the event
		// FIXME pass smeared R9
		int diphoton_id = l.DiphotonCiCSelection(l.phoSUPERTIGHT, l.phoSUPERTIGHT, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
	       
		if (diphoton_id > -1 ) {
		   
		    diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
		    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;
		   
		    float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
		    float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
		    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
		    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
		    TLorentzVector Higgs = lead_p4 + sublead_p4; 	
		   
		    int category = 0; 
		    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
		    if( cur_type != 0 && doMCSmearing ) {
			for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
			    float rewei=1.;
			    float pth = Higgs.Pt();
			    (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), 0. );
			    evweight *= rewei;
			}
		    }
		    float mass = Higgs.M();
			
		    float ptHiggs = Higgs.Pt();
		    // Mass Resolution of the Event
		    massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);

		    double massResolution = massResolutionCalculator->massResolution();

                   // Iterate over each mass point. 
                   for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
                     if (SignalType(cur_type)!=i) continue;

                     // define hypothesis masses for the sidebands
                     float mass_hypothesis = masses[i];
                     float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                     float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                     // define the sidebands
                     float sideband_boundaries[4];
                     sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                     sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                     sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                     sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);
		
                     //Signal Window
                     if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2]){//Signal mass window cut
                       SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolution,mass_hypothesis,evweight);
                       float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada"+names[i] );
                       float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad"+names[i] );
		       categories.push_back(category);
		       bdt_ada_errors.push_back(bdt_ada);
		       bdt_grad_errors.push_back(bdt_grad);
		       weights.push_back(evweight);

		    } else {

		       categories.push_back(-1);
		       bdt_ada_errors.push_back(-100.);
		       bdt_grad_errors.push_back(-100.);
		       weights.push_back(0.);
		    }
		  }

		} else { // In case CiC selection fails now
		   categories.push_back(-1);
		   bdt_ada_errors.push_back(-100.);
		   bdt_grad_errors.push_back(-100.);
		   weights.push_back(0.);
		}

	  }// end loop on systematics steps

	  // Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
          // Iterate over each mass point. 
          for (int i = 2; i<nMasses;i++){ //ignoring masses 105 and 110 for now
              if (SignalType(cur_type)==i){
                  l.rooContainer->InputSystematicSet("sig_BDT_ada"+names[i],(*si)->name(),categories,bdt_ada_errors,weights);
                  l.rooContainer->InputSystematicSet("sig_BDT_grad"+names[i],(*si)->name(),categories,bdt_grad_errors,weights);
              }
	  }
	  // -------------------------------------------------------------------------------------------------------------//
 
	} // Close on Smearers
    } // End Signal Systematic Study 
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
        
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool MvaAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------

/*

double MvaAnalysis::GetPhotonResolution(LoopAll &l,int photon_index){

	// Relevant information
	TLorentzVector *p4    = (TLorentzVector*) (l.pho_p4->At(photon_index));
	TVector3 *sc_xyz = (TVector3*) (l.sc_xyz->At(l.pho_scind[photon_index]));
	double r9 	      = l.pho_r9[photon_index];

  //Make the PhotonFix Classes:
	PhotonFix photonFixer(p4->E(),sc_xyz->Eta(),sc_xyz->Phi(),r9);
	double photonResolution = photonFixer.sigmaEnergy();
	double photonEnergy 	= photonFixer.fixedEnergy();
	// Get the photon-category sigma
	int phoCat = l.PhotonCategory(photon_index,nR9Categories,nEtaCategories);
  	std::string myCategory="";
  	if (eSmearPars.categoryType=="2CatR9_EBEE" )
    	{
      	  if (phoCat==0 || phoCat==1)	myCategory+="EB";
      	  else	myCategory+="EE";
      
          if (phoCat==0 || phoCat==2)   myCategory+="HighR9";
      	  else	myCategory+="LowR9";
    	}
  	else if (eSmearPars.categoryType=="EBEE")
    	{
      	  if (phoCat==0 || phoCat==1)	myCategory+="EB";
      	  else	myCategory+="EE";
    	} 
	// Smearing is applied ON TOP of PhotonFix corrections -> must scale that energy to get addistional resolutions
	double categoryResolution = eSmearPars.smearing_sigma[myCategory]*photonEnergy;	
	return TMath::Sqrt(categoryResolution*categoryResolution + photonResolution*photonResolution);
}

double MvaAnalysis::GetAngleResolutionCorrVtx(LoopAll &l, int lead_index, int sublead_index, int diphoton_index){
  
  TVector3 *vtx_dxdydz = (TVector3*)l.vtx_std_dxdydz->At(l.dipho_vtxind[diphoton_index]);
  double dz = vtx_dxdydz->z();
  
  return PropagateDz(l,lead_index,sublead_index,diphoton_index,dz);
 
}

double MvaAnalysis::GetAngleResolutionWrongVtx(LoopAll &l, int lead_index, int sublead_index, int diphoton_index){
  
  double dz = TMath::Sqrt(2.)*5.8;
  return PropagateDz(l,lead_index,sublead_index,diphoton_index,dz);
}

// propaget vertex resolution to angular resolution
double MvaAnalysis::PropagateDz(LoopAll &l, int lead_index, int sublead_index, int diphoton_index, double dz){

  TLorentzVector *p4_lead = (TLorentzVector*)(l.pho_p4->At(lead_index));
  TLorentzVector *p4_sublead = (TLorentzVector*)(l.pho_p4->At(sublead_index));
  double alpha = p4_lead->Angle(p4_sublead->Vect());
  if (alpha!= p4_sublead->Angle(p4_lead->Vect())) std::cout << "Error: Angle between photons not consistent" << std::endl;
  
  TVector3 *lead_sc_pos = (TVector3*)(l.sc_xyz->At(l.pho_scind[lead_index]));
  TVector3 *sublead_sc_pos = (TVector3*)(l.sc_xyz->At(l.pho_scind[sublead_index]));
  TVector3 *vtx_pos = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_index]);
  
  double x1 = lead_sc_pos->X()-vtx_pos->X();
  double y1 = lead_sc_pos->Y()-vtx_pos->Y();
  double z1 = lead_sc_pos->Z()-vtx_pos->Z();
 
  double x2 = sublead_sc_pos->X()-vtx_pos->X();
  double y2 = sublead_sc_pos->Y()-vtx_pos->Y();
  double z2 = sublead_sc_pos->Z()-vtx_pos->Z();
 
  double r1 = TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  double r2 = TMath::Sqrt(x2*x2+y2*y2+z2*z2);

  double cos_term = TMath::Cos(p4_lead->Eta()-p4_sublead->Eta());
  double sech1 = SecH(p4_lead->Eta());
  double sech2 = SecH(p4_sublead->Eta());
  double tanh1 = TanH(p4_lead->Eta());
  double tanh2 = TanH(p4_sublead->Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2-sech1*sech2*cos_term;

  double ResTerm = (-0.5*dz/denominator)*(numerator1/r1 + numerator2/r2);

  double vertexResolution = ResTerm*2.*(1.-TMath::Cos(alpha))/TMath::Sin(alpha);

  return vertexResolution;

}

// utility funcs for propagating dz
double MvaAnalysis::SecH(double x){
  return 1.0/TMath::CosH(x);
}

double MvaAnalysis::TanH(double x){
  return TMath::TanH(x);
}

*/
int MvaAnalysis::SignalType(int cur_type){
    int i0 = -1;
    if (cur_type == -13 || cur_type == -14 ||  cur_type == -15 || cur_type == -16){//105 
        i0 = 0;}
    else if (cur_type == -17 || cur_type == -18 ||  cur_type == -19 || cur_type == -20){//110 
        i0 = 1;}
    else if (cur_type == -21 || cur_type == -22 ||  cur_type == -23 || cur_type == -24){//115 
        i0 = 2;}
    else if (cur_type == -25 || cur_type == -26 ||  cur_type == -27 || cur_type == -28){//120 
        i0 = 3;}
    else if (cur_type == -37 || cur_type == -38 ||  cur_type == -39 || cur_type == -40){//125 
        i0 = 4;}
    else if (cur_type == -29 || cur_type == -30 ||  cur_type == -31 || cur_type == -32){//130 
        i0 = 5;}
    else if (cur_type == -41 || cur_type == -42 ||  cur_type == -43 || cur_type == -44){//135 
        i0 = 6;}
    else if (cur_type == -33 || cur_type == -34 ||  cur_type == -35 || cur_type == -36){//140 
        i0 = 7;}
    else if (cur_type == -45 || cur_type == -46 ||  cur_type == -47 || cur_type == -48){//145 
        i0 = 8;}
    else if (cur_type == -49 || cur_type == -50 ||  cur_type == -51 || cur_type == -52){//150 
        i0 = 9;}
    else if (cur_type == -53 || cur_type == -54 ||  cur_type == -55 || cur_type == -56){//121 
        i0 = 10;}
    else if (cur_type == -57 || cur_type == -58 ||  cur_type == -59 || cur_type == -60){//123 
        i0 = 11;}
    return i0;
}


void MvaAnalysis::SetBDTInputVariables(TLorentzVector *lead_p4, TLorentzVector *sublead_p4, 
                                       double lead_r9, double sublead_r9, 
                                       double massResolution, double mass_hypothesis, double evweight ){

	TLorentzVector Higgs = *lead_p4 + *sublead_p4; 	
	float mass    = Higgs.M();

    _log_H_pt =  log10( Higgs.Pt());
    _H_eta = fabs(Higgs.Eta());
    _d_phi = fabs(lead_p4->DeltaPhi(*sublead_p4));
    _max_eta = max(fabs(lead_p4->Eta()),fabs(sublead_p4->Eta()));
    _min_r9  = min(lead_r9,sublead_r9);
    _pho1_eta = lead_p4->Eta();
    _pho2_eta = sublead_p4->Eta();

    _mgg = mass;
    _pho1_phi = lead_p4->Phi();
    _pho1_pt = lead_p4->Pt();
    _pho1_r9 = lead_r9;

    _pho2_phi = sublead_p4->Phi();
    _pho2_pt = sublead_p4->Pt();
    _pho2_r9 = sublead_r9;

    _H_pt = Higgs.Pt();
    _Ht = lead_p4->Pt()+sublead_p4->Pt();

    _d_eta = lead_p4->Eta()-sublead_p4->Eta();
    _mod_d_eta = fabs(lead_p4->Eta()-sublead_p4->Eta());
    _cos_theta_star = fabs(lead_p4->E()-sublead_p4->E())/Higgs.P();

    _wt= evweight;

    _pho1_ptOverM = lead_p4->Pt()/mass_hypothesis;
    _pho2_ptOverM = sublead_p4->Pt()/mass_hypothesis;
    _deltaMOverM = (mass-mass_hypothesis)/mass_hypothesis;
    _deltaMOverSigmaM = (mass-mass_hypothesis)/massResolution;
    _sigmaMOverM = massResolution/mass;
    _H_ptOverM    = Higgs.Pt()/mass_hypothesis;

}

