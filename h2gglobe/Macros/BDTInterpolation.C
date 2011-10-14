#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TList.h"
#include "TString.h"
#include "Normalization.C"

using namespace std;

int davidCalls=0;
int fracCalls=0;
int linCalls=0;

TH1F* Interpolate(double massLow, TH1F* low, double massHigh, TH1F* high, double massInt){

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
  name.replace(ind-6,11,strInt);
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

TH1F* linearBin(TH1F *hist){
  assert(hist->GetEntries()!=0);
  assert(hist->GetNbinsX()>2);
  int nBins = hist->GetNbinsX();
  TH1F *temp = new TH1F(Form("lin_%s_%d",hist->GetName(),linCalls),Form("lin_%s_%d",hist->GetName(),linCalls),nBins,0,nBins);
  for (int i=0; i<nBins; i++){
    temp->SetBinContent(i+1,hist->GetBinContent(i+1));
    temp->SetBinError(i+1,hist->GetBinError(i+1));
  }
  linCalls++;
  return temp;
}


void plotDavid(TH1F* bkgT, TH1F* sigT, TH1F* dataT, std::string name){
  
  gROOT->SetBatch();
  system("mkdir -p plots/ada/david");
  system("mkdir -p plots/grad/david");

  std::string bdt;
  TString str = dataT->GetName();
  if (str.Contains("ada")) bdt="ada";
  else if (str.Contains("grad")) bdt="grad";
  else std::cout << "Error find BDT type" << std::endl;
  assert (str.Contains("ada") || str.Contains("grad"));

  if (bkgT->GetNbinsX() != sigT->GetNbinsX() || sigT->GetNbinsX() != dataT->GetNbinsX()) std::cout << "Plot problem: calling plot for histograms with different number of bins" << std::endl;
  assert(bkgT->GetNbinsX() == sigT->GetNbinsX() || sigT->GetNbinsX() == dataT->GetNbinsX());
  
  TH1F *bkg = linearBin(bkgT);
  TH1F *sig = linearBin(sigT);
  TH1F *data = linearBin(dataT);

  gStyle->SetOptStat(0);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();

  TCanvas *canv = new TCanvas();
  bkg->SetLineColor(kBlue);
  bkg->SetFillColor(kBlue-9);
  sig->SetLineColor(kRed);
  sig->Scale(10.);
  data->SetMarkerStyle(20);
  bkg->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  sig->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  data->SetTitle(("Background, data and signal distributions for "+bdt+" "+name).c_str());
  bkg->GetXaxis()->SetTitle("BDT Output Bin");
  sig->GetXaxis()->SetTitle("BDT Output Bin");
  data->GetXaxis()->SetTitle("BDT Output Bin");
  bkg->GetYaxis()->SetTitle("Events");
  sig->GetYaxis()->SetTitle("Events");
  data->GetYaxis()->SetTitle("Events");

  TLegend *leg = new TLegend(0.15,0.6,0.55,0.85);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(bkg,"Background","f");
  leg->AddEntry(data,"Data","lep");
  leg->AddEntry(sig,"Signal (10#timesSM)","l");
  TPaveText *txt = new TPaveText(0.6,0.6,0.8,0.85,"NDC");
  txt->SetFillColor(0);
  txt->SetLineColor(0);
  txt->AddText("#int L = 1.66 fb^{-1}");

  bkg->GetYaxis()->SetRangeUser(0.0,2.*(data->GetMaximum()));

  bkg->Draw("e2");
  sig->Draw("same hist");
  data->Draw("same e");
  leg->Draw("same");
  txt->Draw("same");

  canv->Print(("plots/"+bdt+"/david/"+name+".png").c_str(),"png");

  delete canv;
  delete txt;
  delete leg;
  
  davidCalls++;
}

void plotFrac(TList* HistList, TH1F* compT, std::string name){

  gROOT->SetBatch();
  system("mkdir -p plots/ada/fracs");
  system("mkdir -p plots/grad/fracs");

  std::string bdt;
  TString str = compT->GetName();
  if (str.Contains("ada")) bdt="ada";
  if (str.Contains("grad")) bdt="grad";

  int nHists=HistList->GetEntries();
  TH1F *comp = linearBin(compT);
  TH1F *uandd[nHists];
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TCanvas *canv = new TCanvas("","",700,700);
  
  canv->Divide(1,2,0.1,0.1);
  canv->cd(1);
  gPad->SetPad(0.01,0.3,0.99,0.99);
  gPad->SetBottomMargin(0.0);
  gPad->SetLeftMargin(0.1);
  comp->SetLineColor(1);
  comp->SetFillColor(kGray);
  comp->Draw("hist");
  comp->SetTitle(("Up, down and interpolated templates for "+bdt+" "+name).c_str());
  comp->GetYaxis()->SetTitle("Events");
  comp->GetXaxis()->SetTitle("BDT output bin");
  TLegend *leg = new TLegend(0.2,0.2,0.6,0.8);
  leg->AddEntry(comp,comp->GetName(),"f");
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  
  canv->cd(2);
  gPad->SetPad(0.01,0.01,0.99,0.3);
  gPad->SetTopMargin(0.0);
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.1);
  TF1 *line = new TF1("line","0.0",0.,comp->GetBinLowEdge(comp->GetNbinsX()+1));
  line->SetLineColor(kBlack);
  
  for (int i=0; i<nHists; i++){
    uandd[i] = linearBin((TH1F*)HistList->At(i));
    uandd[i]->SetLineColor((i+1)*2);
    canv->cd(1);
    uandd[i]->DrawCopy("same e");
    uandd[i]->Add(comp,-1);
    uandd[i]->Divide(comp);
    leg->AddEntry(uandd[i],uandd[i]->GetName(),"lep");
    uandd[i]->Scale(100.0);
    canv->cd(2);
    if (i==0) uandd[i]->Draw("e");
    else {
      uandd[i]->Draw("same e");
      line->Draw("same");
    }
    uandd[i]->SetTitle("");
    uandd[i]->GetYaxis()->SetLabelSize(0.08);
    uandd[i]->GetYaxis()->SetRangeUser(-50.,50.);
    uandd[i]->GetYaxis()->SetTitle("#frac{#Delta_{int}}{int} %");
    uandd[i]->GetYaxis()->SetTitleOffset(0.4);
    uandd[i]->GetYaxis()->SetTitleSize(0.08);
    uandd[i]->GetXaxis()->SetLabelSize(0.08);
    uandd[i]->GetXaxis()->SetTitle("BDT output");
    uandd[i]->GetXaxis()->SetTitleSize(0.08);
    canv->cd(1);
   // leg->Draw("same");
  }
  

  canv->Print(("plots/"+bdt+"/fracs/"+name+".png").c_str(),"png");
  delete canv;
  
  fracCalls++;
  
}

std::string getTime(){
  
  std::string timeStr;
  time_t now = time(0);
  struct tm* tm = localtime(&now);
  timeStr = Form("%02d:%02d_%02d-%02d-%02d",tm->tm_hour,tm->tm_min,tm->tm_mday,tm->tm_mon+1,tm->tm_year-100);

  return timeStr;

}

int BDTInterpolation(std::string inFileName,bool Diagnose=false, bool doNorm=true){

  std::cout << getTime() << std::endl;

  system("rm -r plots");
  system("mkdir plots");

  gROOT->SetBatch();
  gStyle->SetOptStat(0);

  // input flags
  bool all=0;
  bool david=0;
  bool frac=0;
 
  if (Diagnose) {
    std::cout << "Diagnostics turned on \n Output plots available in \"plots\" directory \n Diagnostic log available in \"plots/BDTInterpolationDiagnostics.txt\"" << std::endl;
    all = true;
  }
  else std::cout << "Diagnostics turned off" << std::endl;
  if (doNorm) std::cout << "Normalization turned on" << std::endl;
  else std::cout << "Normalization turned off" << std::endl;

  TFile *inFile = new TFile(inFileName.c_str());
  //TFile *inFile = new TFile("/vols/cms02/nw709/hgg/src_cvs/oct13/CMSSW_4_2_8/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/Macros/CMS-HGG_1658pb_mva.root");
  //TFile *inFile = new TFile("RefWorkspaces/CMS-HGG_1658pb_mva.root");
  TFile *outFile = new TFile(Form("%s_interpolated.root",inFileName.c_str()),"RECREATE");
  ofstream diagFile("plots/BDTInterpolationDiagnostics.txt");

  const int nBDTs=2;
  const int nMasses=7;
  std::string BDTtype[nBDTs] = {"ada","grad"};
  std::string BDTmasses[nMasses] = {"115.0","120.0","125.0","130.0","135.0","140.0","150.0"};

// ----- else can just get rebinned histograms straight out of workspace

  std::cout << "Extracting histograms from workspace........." << std::endl;

  diagFile << "Diagnostics for Signal Interpolation run at " << getTime() << std::endl;
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Following orginal histograms rewritten into new workspace: " << std::endl;
  
  // write original histograms in out file
  TList *HistList = inFile->GetListOfKeys();
  for (int j=0; j<HistList->GetSize(); j++){
    TH1F *temp = (TH1F*)inFile->Get(HistList->At(j)->GetName());
    TString name = temp->GetName();
//    if (name.Contains("121") || name.Contains("123") || name.Contains("BDT")) continue;
    std::string tName = temp->GetName();
    for (int i=0; i<nMasses; i++){
      int ind = tName.find(BDTmasses[i]+"_"+BDTmasses[i]);
      if (ind>0) tName.replace(ind,11,BDTmasses[i]);
    }
    temp->SetName(tName.c_str());
      
    outFile->cd();
    temp->Write();
    diagFile << "Histo written: " << temp->GetName() << std::endl;
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
            if (k==-1) orgHistListBelow[bdt][bdtmass]->Add(temp);
            if (k==0) orgHistList[bdt][bdtmass]->Add(temp);
            if (k==1) orgHistListAbove[bdt][bdtmass]->Add(temp);
          }
        }
      }
    }
  }
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Following histo's being used for interpolation: " << std::endl;
  for (int bdt=0; bdt<2; bdt++){
    for (int mass=0; mass<nMasses; mass++){
      diagFile << "BDT: " << BDTtype[bdt] << std::endl;
      diagFile << "Mass: " << BDTmasses[mass] << std::endl;
      for (int syst=0; syst<orgHistList[bdt][mass]->GetSize(); syst++){
        diagFile << "   Central: " << orgHistList[bdt][mass]->At(syst)->GetName() << std::endl;
        if (mass!=0) diagFile << "   Lower:   " << orgHistListBelow[bdt][mass]->At(syst)->GetName() << std::endl;
        if (mass!=nMasses-1) diagFile << "   Upper:   " << orgHistListAbove[bdt][mass]->At(syst)->GetName() << std::endl;
      }
    }
  }

  // now have orgHistList for each BDT and above and below signals.
  // the contain 17 histos - 1 for actual sig and 8 up and down for each systematic

  // ---------- at this points have signal hists for each mass points as well as +- models all rebinned at mass points --------

  std::cout << "Calculating interpolated signal templates........" << std::endl;
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Interpolating intermediate signals from following lower and upper templates" << std::endl;

  TList *orgHistListInt[2][71];
  for (int i=0; i<2; i++) for (int j=0; j<71; j++) orgHistListInt[i][j] = new TList();
  
  int i=0;
  // loop over mass points etc.
  for (double mass=115.; mass<150.5; mass+=0.5){
    if (int(mass*2)%10==0) std::cout << Form("%3.0f",((mass-115.)/35.)*100.) << "% done" << std::endl;
    if (int(mass*2)%10==0 && mass!=145.0) {
      for (int bdt=0; bdt<2; bdt++){  
        TH1F *background = (TH1F*)inFile->Get(Form("th1f_bkg_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        TH1F *data = (TH1F*)inFile->Get(Form("th1f_data_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
        TH1F *signal = (TH1F*)inFile->Get(Form("th1f_sig_%s_%3.1f_%3.1f_cat0",BDTtype[bdt].c_str(),mass,mass));
        std::string name = Form("%3.1f",mass);
        if (all || david) plotDavid(background,signal,data,name);
      }
      continue;
    }
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
      diagFile << "Mass: " << mass << std::endl;
      diagFile << "BDT: " << BDTtype[bdt] << std::endl;
      // loop different histos in list (signal and systematics)
      TH1F *background = (TH1F*)inFile->Get(Form("th1f_bkg_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
      TH1F *data = (TH1F*)inFile->Get(Form("th1f_data_%s_%3.1f_cat0",BDTtype[bdt].c_str(),mass));
      for (int syst=0; syst<orgHistList[bdt][bdtmass]->GetSize(); syst++){
        TH1F *tempSig = (TH1F*)orgHistList[bdt][bdtmass]->At(syst);
        TH1F *tempAbove = (TH1F*)orgHistListAbove[bdt][bdtmass]->At(syst);
        TH1F *tempBelow = (TH1F*)orgHistListBelow[bdt][bdtmass]->At(syst);
        TH1F* tempInt;
        TList* plotList = new TList();
        if (above){
          if (doNorm) tempSig->Scale(1./(GetXsection(double(nearest))*GetBR(double(nearest))));
          if (doNorm) tempAbove->Scale(1./(GetXsection(double(nextNear))*GetBR(double(nextNear))));
          tempInt = Interpolate(double(nearest),tempSig,double(nextNear),tempAbove,mass);
          if (syst==0){
            plotList->Add(tempSig);
            plotList->Add(tempAbove);
            if (all || frac) plotFrac(plotList,tempInt,Form("%3.1f",mass));
          }
          if (doNorm) tempSig->Scale(GetXsection(double(nearest))*GetBR(double(nearest)));
          if (doNorm) tempAbove->Scale(GetXsection(double(nextNear))*GetBR(double(nextNear)));
          double norm = GetNorm(double(nearest),tempSig,double(nextNear),tempAbove,mass);
          diagFile << "   Interpolated:     " << tempInt->GetName() << std::endl;
          diagFile << "   scaled up with:   " << Form("%1.3f",GetXsection(double(nearest))*GetBR(double(nearest))) << std::endl;
          diagFile << "   scaled down with: " << Form("%1.3f",1./norm) << std::endl;
          diagFile << "   scale ratio:      " << Form("%1.3f",(GetXsection(double(nearest))*GetBR(double(nearest)))*1./norm) << std::endl;
          diagFile << "   from:             " << tempSig->GetName() << std::endl;
          diagFile << "   and:              " << tempAbove->GetName() << std::endl;
          if (doNorm) tempInt->Scale(norm/tempInt->Integral());
        }
        else{
          if (doNorm) tempSig->Scale(1./(GetXsection(double(nearest))*GetBR(double(nearest))));
          if (doNorm) tempBelow->Scale(1./(GetXsection(double(nextNear))*GetBR(double(nextNear))));
          tempInt = Interpolate(double(nextNear),tempBelow,double(nearest),tempSig,mass);
          if (syst==0){
            plotList->Add(tempBelow);
            plotList->Add(tempSig);
            if (all || frac) plotFrac(plotList,tempInt,Form("%3.1f",mass));
          }
          if (doNorm) tempSig->Scale(GetXsection(double(nearest))*GetBR(double(nearest)));
          if (doNorm) tempBelow->Scale(GetXsection(double(nextNear))*GetBR(double(nextNear)));
          double norm = GetNorm(double(nextNear),tempBelow,double(nearest),tempSig,mass);
          diagFile << "   Interpolated:     " << tempInt->GetName() << std::endl;
          diagFile << "   scaled up with:   " << Form("%1.3f",GetXsection(double(nearest))*GetBR(double(nearest))) << std::endl;
          diagFile << "   scaled down with: " << Form("%1.3f",1./norm) << std::endl;
          diagFile << "   scale ratio:      " << Form("%1.3f",(GetXsection(double(nearest))*GetBR(double(nearest)))*1./norm) << std::endl;
          diagFile << "   from:             " << tempBelow->GetName() << std::endl;
          diagFile << "   and:              " << tempSig->GetName() << std::endl;
          if (doNorm) tempInt->Scale(norm/tempInt->Integral());
        }
        orgHistListInt[bdt][i]->Add(tempInt);
        std::string name = Form("%3.1f",mass);
        if (syst==0 && (all || david)) plotDavid(background,tempInt,data,name);
        delete plotList;
      }
    }
    i++;
  }
 
  outFile->cd();
  
  diagFile << "---------------------------------------------------------------------" << std::endl;
  diagFile << "Writing following interpolated histograms to file: " << std::endl;
  for (int l=0; l<2; l++) for (int j=0; j<71; j++) for (int k=0; k<orgHistListInt[l][j]->GetSize(); k++) {
    orgHistListInt[l][j]->At(k)->Write();
    diagFile << orgHistListInt[l][j]->At(k)->GetName() << std::endl;
  }
  
  TList* endList = outFile->GetListOfKeys();
  for (double mass=115.0; mass<150.5; mass+=0.5){
    diagFile << mass << std::endl;
    std::pair<int,int> nearestPair = findNearest(mass);
    double nearest = nearestPair.first;
    double nextNear = nearestPair.second;
    for (int k=0; k<endList->GetSize(); k++){
      TString hName = endList->At(k)->GetName();
      if (hName.Contains("sigma") || !hName.Contains("sig")) continue;
      if (hName.Contains(Form("a_%3.1f_cat0",nearest)) || hName.Contains(Form("d_%3.1f_cat0",nearest)) || hName.Contains(Form("%3.1f_%3.1f_cat0",nearest,nextNear)) || hName.Contains(Form("%3.1f_cat0",mass))) {
        TH1F *temp = (TH1F*)outFile->Get(hName);
        diagFile << "    " << setw(30) << temp->GetName();
        diagFile << "    " << setw(8) << temp->GetEntries();
        diagFile << "    " << setw(8) << temp->Integral() << std::endl;
      }
    }
  }
    
  outFile->Close();

  std::cout << "Diagnostics log written to \"BDTInterpolationDiagnostics.txt\"" << std::endl;
  if (davidCalls>0) std::cout << davidCalls << " plots written to: plots/ada/david/ \n                  plots/grad/david/" << std::endl;
  if (fracCalls>0) std::cout << fracCalls << " plots written to: plots/ada/frac/ \n                  plots/grad/frac/" << std::endl;

  std::string TIME_DATE = getTime();
  if (Diagnose){
    system(("python make_html.py "+TIME_DATE).c_str());
    system(("mkdir -p ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str());
    system(("cp -r plots ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str());
    system("whoami > temp.txt");
    std::ifstream temp("temp.txt");
    std::string result;
    temp >> result;
    temp.close();
    system("rm temp.txt");
    std::cout << ("Plots avaiable to view in ~/public_html/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/").c_str() << std::endl;
    std::cout << "If working on /vols/ at IC plots avaliable to view at www.hep.ph.ic.ac.uk/~"+result+"/h2g/MVA/SigInt/Diagnostics/"+TIME_DATE+"/plots/plots.html" << std::endl;
  }
  std::cout << "New workspace written to " << outFile->GetName() << std::endl;
  return 0;
}

