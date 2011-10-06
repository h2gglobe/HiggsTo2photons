#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TList.h"
#include "TString.h"
#include "Normalization.C"

std::vector<double> optimizedBinning(TH1F *hb, int nTargetBins, bool revise_target, bool use_n_target){
	// Return a set of bins which are "smoother" 

	if (revise_target) {
		if (use_n_target){
		   std::cerr << "WARNING -- RooContainer::OptimizedBinning -- Can't use number of Entries as target in revised binning algo " << std::endl; 
		   use_n_target = false;  // geometric algo always use revised number of bins, not number of entries
		
		}
	}

	int nBins = hb->GetNbinsX();
	std::vector<double> binEdges;
	binEdges.push_back(hb->GetBinLowEdge(1));

	double targetNumbers;
	if (use_n_target) targetNumbers = nTargetBins; 
	else targetNumbers = hb->Integral()/nTargetBins;

	if (hb->Integral() < 2*targetNumbers){
		std::cout << "RooContainer::OptimizedBinning -- Not enough entries in histogram for target numbers calculated - " 
			  << targetNumbers 
			  << ", Returning current bin boundaries "  << std::endl;
		for (int j=2;j<=nBins+1;j++) binEdges.push_back(hb->GetBinLowEdge(j));
		return binEdges;
	}

	double sumBin = 0;
	int i=1;
	while (i<=nBins){
		if (revise_target) targetNumbers = hb->Integral(i,nBins)/nTargetBins;
		sumBin=hb->GetBinContent(i);
		double highEdge=hb->GetBinLowEdge(i+1);

		bool carryOn = sumBin <= targetNumbers;
		while ( carryOn){
			if (i<nBins){
			  sumBin+=hb->GetBinContent(i+1);
			  highEdge = hb->GetBinLowEdge(i+2);
			  carryOn =(sumBin <targetNumbers && i<=nBins);
			  i++;
			} else {
			  highEdge = hb->GetBinLowEdge(i+1);
			  carryOn=false;
			}
		}
	        binEdges.push_back(highEdge);
		i++;
	}
        if (sumBin < 10) binEdges.erase(binEdges.end()-2);
	return binEdges;

}
// ----------------------------------------------

TH1F* rebinBinnedDataset(std::string new_name,std::string name, TH1F *hb,std::vector<double> binEdges){

	double *arrBins = new double[binEdges.size()];
	int j=0;
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		arrBins[j]=*it;
		j++;
		
	}
	//const char *h_name = (const char *) hb->GetName;
	//const char *title  = (const char *) hb->GetTitle;
	
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,hb->GetName(),arrBins);
	hbnew->SetName(Form("th1f_%s",new_name.c_str()));
	//cout << "title for new re-binned histogream - " << hb->GetTitle()<<endl; 
	hbnew->SetTitle(hb->GetTitle());

  return hbnew;
}
// ----------------------------------------------

TH1F* Interpolate(TH1F* low, TH1F* high, double massLow, double massHigh, double massInt){

  if (low->GetNbinsX()!=high->GetNbinsX()) std::cout << "Cannot interpolate differently binned histograms" << std::endl;
  assert(low->GetNbinsX()==high->GetNbinsX());

  TH1F *interp = (TH1F*)low->Clone();
  for (int i=0; i<interp->GetNbinsX(); i++){
    double OutLow = low->GetBinContent(i+1);
    double OutHigh = high->GetBinContent(i+1);
    double OutInt = (OutHigh*(massInt-massLow)-OutLow*(massInt-massHigh))/(massHigh-massLow);
    interp->SetBinContent(i+1,OutInt);
  }

  std::string name = low->GetName();
  std::string strLow = Form("%3.1f",massLow);
  std::string strInt = Form("%3.1f",massInt);
  int ind = name.rfind(strLow);
  name.replace(ind,5,strInt);
  interp->SetName(name.c_str());

  return interp;
   
}

std::pair<int,int> findNearest(double mass){
  
  double bestdiff=1000.0;
  double nbestdiff=1000.0;
  int nearest;
  int nnearest;

  for (int i=115; i<155; i+=5){
    if (i==145) continue;
    double diff = std::fabs((mass+0.001)-double(i));
    if (diff<bestdiff) {
      bestdiff=diff;
      nearest = i;
    }
  }
  for (int i=115; i<155; i+=5){
    if (i==145) continue;
    double diff = std::fabs((mass+0.001)-double(i));
    if (diff>bestdiff && diff<nbestdiff){
      nbestdiff=diff;
      nnearest=i;
    }
  }
  std::pair<int,int> nearestPair(nearest,nnearest);
  return nearestPair;
}

int getIndex(int mass){

  int index;
  if (mass==115) index=0;
  if (mass==120) index=1;
  if (mass==125) index=2;
  if (mass==130) index=3;
  if (mass==135) index=4;
  if (mass==140) index=5;
  if (mass==150) index=6;

  return index;
}

int main(){
  
  TFile *inFile = new TFile("/vols/cms02/h2g/interp_workspace/CMS-HGG_1658pb_mva.root");
  TFile *outFile = new TFile("/vols/cms02/h2g/interp_workspace/CMS-HGG_1658pb_mva_int_02-10-11.root","RECREATE");

  const int nBDTs=2;
  const int nMasses=7;
  std::string BDTtype[nBDTs] = {"ada","grad"};
  std::string BDTmasses[nMasses] = {"115.0","120.0","125.0","130.0","135.0","140.0","150.0"};

  std::vector<double> optBins[2][7];

  TH1F *sigHist[2][7];
  TH1F *sigHistAbove[2][7];
  TH1F *sigHistBelow[2][7];

// ---- if rebinning required can use this step -------
/*
  // get 14 optimized binnings for each mass point
  for (int bdt=0; bdt<nBDTs; bdt++){
    for (int bdtmass=0; bdtmass<nMasses; bdtmass++){
      TH1F *bkg = (TH1F*)inFile->Get(("th1f_bkg_BDT_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_cat0").c_str());
      optBins[bdt][bdtmass] = optimizedBinning(bkg,50,false,true);
      TH1F *sig = (TH1F*)inFile->Get(("th1f_sig_BDT_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_cat0").c_str());
      // rebin signal
      sigHist[bdt][bdtmass] = rebinBinnedDataSet(("sig_reb_"+BDTmasses[bdtmass]+"_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]).c_str()," ", sig, optBins[bdt][bdtmass]);
      // rebin above and below
      if (i>0) {
        TH1F *below = (TH1F*)inFile->Get(("th1f_sig_BDT_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_"+BDTmasses[bdtmass-1]+"_cat0").c_str());
        sigHistBelow[bdt][bdtmass] = rebinBinnedDataSet(("sig_bel_reb_"+BDTmasses[bdtmass]+"_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]).c_str()," ", below, optBins[bdt][bdtmass]);
      }
      if (i<nMasses-1) {
        TH1F *above = (TH1F*)inFile->Get(("th1f_sig_BDT_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_"+BDTmasses[bdtmass+1]+"_cat0").c_str());
        sigHistAbove[bdt][bdtmass] = rebinBinnedDataSet(("sig_abo_reb_"+BDTmasses[bdtmass]+"_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]).c_str()," ", below, optBins[bdt][bdtmass]);
      }
    }
  }
*/

// ----- else can just get rebinned histograms straight out of workspace

  // write original histograms in out file
  TList *HistList = inFile->GetListOfKeys();
  outFile->cd();
  for (int j=0; j<HistList->GetSize(); j++){
    TH1F *temp = (TH1F*)inFile->Get(HistList->At(j)->GetName());
    temp->Write();
  }
  // get lists of middle, upper and lower templates for each mass
  TList *orgHistList[2][7];
  TList *orgHistListBelow[2][7];
  TList *orgHistListAbove[2][7];
  for (int i=0; i<2; i++) {
    for (int j=0; j<7; j++) {
      orgHistList[i][j]=new TList();
      orgHistListBelow[i][j]=new TList();
      orgHistListAbove[i][j]=new TList();
    }
  }
  for (int j=0; j<HistList->GetSize(); j++){
    TString HistName(HistList->At(j)->GetName());
    for (int bdt=0; bdt<nBDTs; bdt++){
      for (int bdtmass=0; bdtmass<nMasses; bdtmass++){
        for (int k=-1; k<2; k++){
          if ((bdtmass==0 && k==-1) || (bdtmass==nMasses-1 && k==1)) continue;
          if (HistName.Contains(("sig_"+BDTtype[bdt]+"_"+BDTmasses[bdtmass]+"_"+BDTmasses[bdtmass+k]).c_str())){
            TH1F *temp = (TH1F*)inFile->Get(HistName.Data());
            // scale by 1/XS*BR
            double XS = GetXsection(std::atof(BDTmasses[bdtmass+k].c_str()));
            double BR = GetBR(std::atof(BDTmasses[bdtmass+k].c_str()));
            temp->Scale(1./(XS*BR));
            if (k==-1) orgHistListBelow[bdt][bdtmass]->Add(temp);
            if (k==0) orgHistList[bdt][bdtmass]->Add(temp);
            if (k==1) orgHistListAbove[bdt][bdtmass]->Add(temp);

          }
        }
      }
    }
  }
  // now have orgHistList for each BDT and above and below signals.
  // the contain 17 histos - 1 for actual sig and 8 up and down for each systematic

  // ---------- at this points have signal hists for each mass points as well as +- models all rebinned at mass points --------

  TList *orgHistListInt[2][71];
  for (int i=0; i<2; i++) for (int j=0; j<71; j++) orgHistListInt[i][j] = new TList();
  
  int i=0;
  // loop over mass points etc.
  for (double mass=115.; mass<150.5; mass+=0.5){
    //std::cout << i << std::endl; 
    if (int(mass*2)%10==0 && mass!=145.0) continue;
    std::pair<int,int> nearestPair = findNearest(mass);
    bool above;
    int nearest = nearestPair.first;
    int nextNear =nearestPair.second;
    if (nearest-nextNear < 0) above = true;
    else above = false;
    
    //std::cout << mass << " bracketed by: " << nearestPair.first << " " << nearestPair.second << " " << above << std::endl;

    int bdtmass = getIndex(nearest); // gives index of bdt to use

    // loop bdt type
    for (int bdt=0; bdt<2; bdt++){
      // loop different histos in list (signal and systematics)
      for (int syst=0; syst<orgHistList[bdt][bdtmass]->GetSize(); syst++){
        TH1F *tempSig = (TH1F*)orgHistList[bdt][bdtmass]->At(syst);
        TH1F *tempAbove = (TH1F*)orgHistListAbove[bdt][bdtmass]->At(syst);
        TH1F *tempBelow = (TH1F*)orgHistListBelow[bdt][bdtmass]->At(syst);
        TH1F* tempInt;
        if (above){
          tempInt = Interpolate(tempSig,tempAbove, double(nearest),double(nextNear),mass);
          double norm = GetNorm(double(nearest),tempSig,double(nextNear),tempAbove,mass);
          tempInt->Scale(norm);
        }
        else{
          tempInt = Interpolate(tempBelow,tempSig, double(nextNear),double(nearest),mass);
          double norm = GetNorm(double(nextNear),tempBelow,double(nearest),tempSig,mass);
          tempInt->Scale(norm);
        }
        orgHistListInt[bdt][i]->Add(tempInt);
      }
    }
    i++;
  }
 
  outFile->cd();
  
  for (int i=0; i<2; i++) for (int j=0; j<71; j++) for (int k=0; k<orgHistListInt[i][j]->GetSize(); k++) {
    orgHistListInt[i][j]->At(k)->Write();
  }

  return 0;
}

