#include "TMath.h"
#include <iomanip>

// standard includes
#include <cmath>
#include <map>
#include <set>
#include <vector>

void backgroundModel_shapeBias(int mass_in=120, bool www=false, TString outdirname="BDTplots_all", bool zoom=false) {

  bool rebin=true;
  bool scale_nentries=true;

  float sidebandWidth=0.02;
  TString sb_str="_2pc";
  TString sob_str="_sob";

  float signalRegionWidth=0.07;
  const float dataMC_sf = 1.;

  TString mass_str;
  mass_str+=mass_in;
  TString mass_str2= mass_str;
  if (mass_in==150) mass_str2="145";

  float mass_sb[7];

  mass_sb[2] = float(mass_in)*(1-signalRegionWidth)/(1+sidebandWidth);
  mass_sb[1] = mass_sb[2]*(1-sidebandWidth)/(1+sidebandWidth);
  mass_sb[0] = mass_sb[1]*(1-sidebandWidth)/(1+sidebandWidth);
  mass_sb[3] = float(mass_in);
  mass_sb[4] = float(mass_in)*(1+signalRegionWidth)/(1-sidebandWidth);
  mass_sb[5] = mass_sb[4]*(1+sidebandWidth)/(1-sidebandWidth);
  mass_sb[6] = mass_sb[5]*(1+sidebandWidth)/(1-sidebandWidth);

  TString outdir;
  if (www) {
    outdir = "/afs/cern.ch/user/f/futyand/www/hgg/"+outdirname+"/";
  } else {
    outdir = outdirname+"/";
  }

  int nbins[2];
  Double_t xbins_nominal_grad[20];
  Double_t xbins_nominal_ada[20];

  TH1* hist_nominal_grad;
  TH1* hist_nominal_ada;

  //TFile *f_bdtout_nominal = TFile::Open("CMS-HGG_4721pb"+sb_str+sob_str+".root");
  //TFile *f_bdtout_nominal = TFile::Open("/tmp/futyand/CMS-HGG_4700pb_02-12-11.root");
  TFile *f_bdtout_nominal = TFile::Open("CMS-HGG_4721pb_9Dec_binning.root");
  hist_nominal_grad = (TH1*)(f_bdtout_nominal->Get("th1f_data_grad_"+mass_str+".0_cat0"))->Clone();
  hist_nominal_ada = (TH1*)(f_bdtout_nominal->Get("th1f_data_ada_"+mass_str+".0_cat0"))->Clone();

  nbins[0] = hist_nominal_grad->GetNbinsX();
  for (int ibin=0; ibin<nbins[0]+1; ibin++) {
    xbins_nominal_grad[ibin] = hist_nominal_grad->GetBinLowEdge(ibin+1);
  }
  nbins[1] = hist_nominal_ada->GetNbinsX();
  for (int ibin=0; ibin<nbins[1]+1; ibin++) {
    xbins_nominal_ada[ibin] = hist_nominal_ada->GetBinLowEdge(ibin+1);
  }

  f_bdtout_nominal->Close();

  TH1* hist_data[15][2][3];
  TH1* hist_data_ref[2][3];
  TH1* hist_bkg[15][2][3];
  TH1* hist_data_fine[15][2][3];
  TH1* hist_data_fine_ref[2][3];
  TH1* hist_bkg_fine[15][2][3];
  TH1* hist_data_sidebands[7][3];
  float bias_sf[3];

  TFile *f_bdtout[3];

  for (int loose=0; loose<3; loose++) {

    TString loose_str;
    if (loose==0) {
      loose_str="";
    } else if (loose==1) {
      loose_str="_loose";
    } else if (loose==2) {
      loose_str="_veryloose";
    } else {
      cout << "loose must be 0, 1 or 2" << endl;
    }

    //f_bdtout[loose] = TFile::Open("CMS-HGG_4721pb"+sb_str+loose_str+sob_str+".root");
    f_bdtout[loose] = TFile::Open("CMS-HGG_4721pb_9Dec"+loose_str+".root");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    gStyle->SetCanvasColor(0);
    gStyle->SetFrameBorderMode(0);

    f_bdtout[loose]->cd();

    hist_data_fine[0][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3low_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[1][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2low_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[2][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1low_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[3][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3low_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[4][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2low_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[5][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1low_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[6][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[7][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[8][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_grad_115.0_cat0"))->Clone();
    hist_data_fine[9][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[10][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[11][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_grad_130.0_cat0"))->Clone();
    hist_data_fine[12][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_grad_150.0_cat0"))->Clone();
    hist_data_fine[13][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_grad_150.0_cat0"))->Clone();
    hist_data_fine[14][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_grad_150.0_cat0"))->Clone();
    hist_data_fine_ref[0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_BDT_grad_"+mass_str+".0_cat0"))->Clone();

    hist_bkg_fine[0][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3low_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[1][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2low_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[2][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1low_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[3][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3low_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[4][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2low_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[5][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1low_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[6][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[7][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[8][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_grad_115.0_cat0"))->Clone();
    hist_bkg_fine[9][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[10][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[11][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_grad_130.0_cat0"))->Clone();
    hist_bkg_fine[12][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_grad_150.0_cat0"))->Clone();
    hist_bkg_fine[13][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_grad_150.0_cat0"))->Clone();
    hist_bkg_fine[14][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_grad_150.0_cat0"))->Clone();

    hist_data_fine[0][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3low_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[1][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2low_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[2][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1low_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[3][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3low_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[4][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2low_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[5][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1low_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[6][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[7][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[8][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_ada_115.0_cat0"))->Clone();
    hist_data_fine[9][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[10][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[11][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_ada_130.0_cat0"))->Clone();
    hist_data_fine[12][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_1high_BDT_ada_150.0_cat0"))->Clone();
    hist_data_fine[13][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_2high_BDT_ada_150.0_cat0"))->Clone();
    hist_data_fine[14][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_3high_BDT_ada_150.0_cat0"))->Clone();
    hist_data_fine_ref[1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_BDT_ada_"+mass_str+".0_cat0"))->Clone();

    hist_bkg_fine[0][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3low_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[1][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2low_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[2][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1low_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[3][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3low_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[4][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2low_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[5][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1low_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[6][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[7][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[8][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_ada_115.0_cat0"))->Clone();
    hist_bkg_fine[9][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[10][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[11][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_ada_130.0_cat0"))->Clone();
    hist_bkg_fine[12][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_BDT_ada_150.0_cat0"))->Clone();
    hist_bkg_fine[13][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_BDT_ada_150.0_cat0"))->Clone();
    hist_bkg_fine[14][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_BDT_ada_150.0_cat0"))->Clone();

    hist_data_fine[0][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_3low_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[1][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_2low_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[2][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_1low_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[3][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_3low_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[4][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_2low_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[5][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_1low_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[6][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_1high_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[7][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_2high_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[8][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_3high_grad_115.0_cat0",xbins_nominal_grad);
    hist_data_fine[9][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_1high_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[10][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_2high_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[11][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_3high_grad_130.0_cat0",xbins_nominal_grad);
    hist_data_fine[12][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_1high_grad_150.0_cat0",xbins_nominal_grad);
    hist_data_fine[13][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_2high_grad_150.0_cat0",xbins_nominal_grad);
    hist_data_fine[14][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_3high_grad_150.0_cat0",xbins_nominal_grad);
    hist_data_fine_ref[0][loose]->Rebin(nbins[0],"th1f_nominal_data_grad_"+mass_str+".0_cat0",xbins_nominal_grad);

    hist_bkg_fine[0][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_3low_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[1][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_2low_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[2][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_1low_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[3][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_3low_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[4][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_2low_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[5][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_1low_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[6][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_1high_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[7][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_2high_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[8][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_3high_grad_115.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[9][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_1high_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[10][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_2high_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[11][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_3high_grad_130.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[12][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_1high_grad_150.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[13][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_2high_grad_150.0_cat0",xbins_nominal_grad);
    hist_bkg_fine[14][0][loose]->Rebin(nbins[0],"th1f_nominal_bkg_mc_3high_grad_150.0_cat0",xbins_nominal_grad);

    hist_data_fine[0][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_3low_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[1][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_2low_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[2][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_1low_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[3][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_3low_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[4][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_2low_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[5][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_1low_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[6][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_1high_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[7][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_2high_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[8][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_3high_ada_115.0_cat0",xbins_nominal_ada);
    hist_data_fine[9][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_1high_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[10][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_2high_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[11][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_3high_ada_130.0_cat0",xbins_nominal_ada);
    hist_data_fine[12][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_1high_ada_150.0_cat0",xbins_nominal_ada);
    hist_data_fine[13][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_2high_ada_150.0_cat0",xbins_nominal_ada);
    hist_data_fine[14][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_3high_ada_150.0_cat0",xbins_nominal_ada);
    hist_data_fine_ref[1][loose]->Rebin(nbins[1],"th1f_nominal_data_ada_"+mass_str+".0_cat0",xbins_nominal_ada);

    hist_bkg_fine[0][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_3low_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[1][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_2low_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[2][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_1low_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[3][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_3low_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[4][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_2low_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[5][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_1low_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[6][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_1high_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[7][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_2high_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[8][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_3high_ada_115.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[9][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_1high_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[10][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_2high_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[11][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_3high_ada_130.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[12][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_1high_ada_150.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[13][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_2high_ada_150.0_cat0",xbins_nominal_ada);
    hist_bkg_fine[14][1][loose]->Rebin(nbins[1],"th1f_nominal_bkg_mc_3high_ada_150.0_cat0",xbins_nominal_ada);

    hist_data[0][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3low_grad_115.0_cat0"))->Clone();
    hist_data[1][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2low_grad_115.0_cat0"))->Clone();
    hist_data[2][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1low_grad_115.0_cat0"))->Clone();
    hist_data[3][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3low_grad_130.0_cat0"))->Clone();
    hist_data[4][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2low_grad_130.0_cat0"))->Clone();
    hist_data[5][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1low_grad_130.0_cat0"))->Clone();
    hist_data[6][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_grad_115.0_cat0"))->Clone();
    hist_data[7][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_grad_115.0_cat0"))->Clone();
    hist_data[8][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_grad_115.0_cat0"))->Clone();
    hist_data[9][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_grad_130.0_cat0"))->Clone();
    hist_data[10][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_grad_130.0_cat0"))->Clone();
    hist_data[11][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_grad_130.0_cat0"))->Clone();
    hist_data[12][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_grad_150.0_cat0"))->Clone();
    hist_data[13][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_grad_150.0_cat0"))->Clone();
    hist_data[14][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_grad_150.0_cat0"))->Clone();
    hist_data_ref[0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_data_grad_"+mass_str+".0_cat0"))->Clone();

    hist_bkg[0][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3low_grad_115.0_cat0"))->Clone();
    hist_bkg[1][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2low_grad_115.0_cat0"))->Clone();
    hist_bkg[2][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1low_grad_115.0_cat0"))->Clone();
    hist_bkg[3][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3low_grad_130.0_cat0"))->Clone();
    hist_bkg[4][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2low_grad_130.0_cat0"))->Clone();
    hist_bkg[5][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1low_grad_130.0_cat0"))->Clone();
    hist_bkg[6][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_grad_115.0_cat0"))->Clone();
    hist_bkg[7][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_grad_115.0_cat0"))->Clone();
    hist_bkg[8][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_grad_115.0_cat0"))->Clone();
    hist_bkg[9][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_grad_130.0_cat0"))->Clone();
    hist_bkg[10][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_grad_130.0_cat0"))->Clone();
    hist_bkg[11][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_grad_130.0_cat0"))->Clone();
    hist_bkg[12][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_grad_150.0_cat0"))->Clone();
    hist_bkg[13][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_grad_150.0_cat0"))->Clone();
    hist_bkg[14][0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_grad_150.0_cat0"))->Clone();

    hist_data[0][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3low_ada_115.0_cat0"))->Clone();
    hist_data[1][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2low_ada_115.0_cat0"))->Clone();
    hist_data[2][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1low_ada_115.0_cat0"))->Clone();
    hist_data[3][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3low_ada_130.0_cat0"))->Clone();
    hist_data[4][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2low_ada_130.0_cat0"))->Clone();
    hist_data[5][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1low_ada_130.0_cat0"))->Clone();
    hist_data[6][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_ada_115.0_cat0"))->Clone();
    hist_data[7][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_ada_115.0_cat0"))->Clone();
    hist_data[8][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_ada_115.0_cat0"))->Clone();
    hist_data[9][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_ada_130.0_cat0"))->Clone();
    hist_data[10][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_ada_130.0_cat0"))->Clone();
    hist_data[11][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_ada_130.0_cat0"))->Clone();
    hist_data[12][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_1high_ada_150.0_cat0"))->Clone();
    hist_data[13][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_2high_ada_150.0_cat0"))->Clone();
    hist_data[14][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_3high_ada_150.0_cat0"))->Clone();
    hist_data_ref[1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_data_ada_"+mass_str+".0_cat0"))->Clone();

    hist_bkg[0][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3low_ada_115.0_cat0"))->Clone();
    hist_bkg[1][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2low_ada_115.0_cat0"))->Clone();
    hist_bkg[2][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1low_ada_115.0_cat0"))->Clone();
    hist_bkg[3][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3low_ada_130.0_cat0"))->Clone();
    hist_bkg[4][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2low_ada_130.0_cat0"))->Clone();
    hist_bkg[5][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1low_ada_130.0_cat0"))->Clone();
    hist_bkg[6][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_ada_115.0_cat0"))->Clone();
    hist_bkg[7][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_ada_115.0_cat0"))->Clone();
    hist_bkg[8][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_ada_115.0_cat0"))->Clone();
    hist_bkg[9][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_ada_130.0_cat0"))->Clone();
    hist_bkg[10][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_ada_130.0_cat0"))->Clone();
    hist_bkg[11][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_ada_130.0_cat0"))->Clone();
    hist_bkg[12][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_1high_ada_150.0_cat0"))->Clone();
    hist_bkg[13][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_2high_ada_150.0_cat0"))->Clone();
    hist_bkg[14][1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_nominal_bkg_mc_3high_ada_150.0_cat0"))->Clone();

    hist_data_sidebands[0][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[1][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[2][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[3][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_data_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[4][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_1high_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[5][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_2high_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data_sidebands[6][loose] = (TH1*)(f_bdtout[loose]->Get("th1f_bkg_3high_grad_"+mass_str+".0_cat0"))->Clone();

    float nentries_sb[7];
    for (int i=0; i<7; i++) nentries_sb[i] = hist_data_sidebands[i][loose]->Integral();
  
    bias_sf[loose]=0.;
    int sb;
    for (int i=0; i<7; i++) {
      if (i!=3) {
	float bias_sfi = nentries_sb[i]/nentries_sb[3]*(mass_sb[i]-mass_sb[3]);
	bias_sf[loose]+=bias_sfi;
      }
    }
    bias_sf[loose]/=6;
    cout << "bias_sf[" << loose << "] = " << bias_sf[loose] << endl;
  }

  TCanvas *canvas[2];
  TString gradada_str[2]={"_grad","_ada"};

  TGraph *G_bdtout_data[10][2][3];
  TGraph *G_bdtout_mc[10][2][3];
  TGraph *G_gradient_data[2][3];
  TGraph *G_gradient_mc[2][3];
  TGraph *G_bias_data[2][3];
  TGraph *G_bias_mc[2][3];
  TMultiGraph *mG_gradient_data[2];
  TMultiGraph *mG_gradient_mc[2];
  TMultiGraph *mG_gradient[2];
  TMultiGraph *mG_gradient_loose[2];
  TMultiGraph *mG_gradient_veryloose[2];
  TMultiGraph *mG_bias[2];
  TMultiGraph *mG_bias_loose[2];
  TMultiGraph *mG_bias_veryloose[2];
  float gradient_data[10][2][3];
  float gradient_mc[10][2][3];
  float gradient_data_err[10][2][3];
  float gradient_mc_err[10][2][3];
  float gradient_data_temp[10];
  float gradient_mc_temp[10];
  float gradient_data_err_temp[10];
  float gradient_mc_err_temp[10];
  float bias_data[10];
  float bias_mc[10];
  float bias_data_err[10];
  float bias_mc_err[10];

  float mass[15];
  float mass_err[15];
  for(int i=0; i<15; i++) mass_err[i]=0.;
  float bdtout_data[15];
  float bdtout_data_err[15];
  float bdtout_mc[15];
  float bdtout_mc_err[15];
  float bin[10];
  float bin_err[10];
  for(int i=0; i<10; i++) {
    bin[i]=0.5+float(i);
    bin_err[i]=0.;
  }
  float bin_shifted[10];
  float binshift[3]={-0.2,0.,0.2};

  mass[2] = 115.*(1-signalRegionWidth)/(1+sidebandWidth);
  mass[1] = mass[2]*(1-sidebandWidth)/(1+sidebandWidth);
  mass[0] = mass[1]*(1-sidebandWidth)/(1+sidebandWidth);
  mass[5] = 130.*(1-signalRegionWidth)/(1+sidebandWidth);
  mass[4] = mass[5]*(1-sidebandWidth)/(1+sidebandWidth);
  mass[3] = mass[4]*(1-sidebandWidth)/(1+sidebandWidth);
  mass[6] = 115.*(1+signalRegionWidth)/(1-sidebandWidth);
  mass[7] = mass[6]*(1+sidebandWidth)/(1-sidebandWidth);
  mass[8] = mass[7]*(1+sidebandWidth)/(1-sidebandWidth);
  mass[9] = 130.*(1+signalRegionWidth)/(1-sidebandWidth);
  mass[10] = mass[9]*(1+sidebandWidth)/(1-sidebandWidth);
  mass[11] = mass[10]*(1+sidebandWidth)/(1-sidebandWidth);
  mass[12] = 150.*(1+signalRegionWidth)/(1-sidebandWidth);
  mass[13] = mass[12]*(1+sidebandWidth)/(1-sidebandWidth);
  mass[14] = mass[13]*(1+sidebandWidth)/(1-sidebandWidth);

  float sideband_boundaries[2][15];
  for(int i=0; i<15; i++) {
    //cout << mass[i] << endl;
    sideband_boundaries[0][i] = mass[i]*(1-signalRegionWidth);
    sideband_boundaries[1][i] = mass[i]*(1+signalRegionWidth);
  }

  Double_t nentries[2];
  TH1* hist_bkg_scaled[15][2][3];
  TH1* hist_data_scaled[15][2][3];

  for (int loose=0; loose<3; loose++) {
    for(int j=0; j<2; j++) {
      if (!scale_nentries) hist_data_ref[j][loose]->Scale(1./hist_data_ref[j][loose]->Integral());
      nentries[j] = hist_data_ref[j][loose]->Integral();
      for(int i=0; i<15; i++) {
	hist_data_scaled[i][j][loose]=(TH1*)hist_data[i][j][loose]->Clone();
	hist_data_scaled[i][j][loose]->Scale(nentries[j]/hist_data[i][j][loose]->Integral());
	hist_bkg_scaled[i][j][loose]=(TH1*)hist_bkg[i][j][loose]->Clone();
	hist_bkg_scaled[i][j][loose]->Scale(nentries[j]/hist_bkg[i][j][loose]->Integral());
      }
    }
  }

  TLine *line[2];

  for(int j=0; j<2; j++) {

    line[j] = new TLine(0.,0.,float(nbins[j]),0.);
    line[j]->SetLineColor(4);
    line[j]->SetLineWidth(2);

    canvas[j] = new TCanvas("canvas"+gradada_str[j],"canvas"+gradada_str[j],1600,1000);
    canvas[j]->Divide(nbins[j],6);
    canvas[j]->SetFillColor(0);

    for (int loose=0; loose<3; loose++) {
      for(int bdtbin=0; bdtbin<nbins[j]; bdtbin++) {
	for(int i=0; i<15; i++) {
	  bdtout_data[i] = hist_data_scaled[i][j][loose]->GetBinContent(bdtbin+1);
	  bdtout_data_err[i] = hist_data_scaled[i][j][loose]->GetBinError(bdtbin+1);
	  bdtout_mc[i] = hist_bkg_scaled[i][j][loose]->GetBinContent(bdtbin+1);
	  bdtout_mc_err[i] = hist_bkg_scaled[i][j][loose]->GetBinError(bdtbin+1);
	}

	canvas[j]->cd(bdtbin+1+(nbins[j]*loose));
	G_bdtout_data[bdtbin][j][loose] = new TGraphErrors(15,mass,bdtout_data,mass_err,bdtout_data_err);
	G_bdtout_data[bdtbin][j][loose]->SetMarkerSize(.8);
	G_bdtout_data[bdtbin][j][loose]->SetMarkerColor(1);
	G_bdtout_data[bdtbin][j][loose]->SetMarkerStyle(20+loose);
	G_bdtout_data[bdtbin][j][loose]->Draw("AP");
	G_bdtout_data[bdtbin][j][loose]->GetXaxis()->SetLabelSize(0.07);
	G_bdtout_data[bdtbin][j][loose]->GetYaxis()->SetLabelSize(0.07);
	G_bdtout_data[bdtbin][j][loose]->GetXaxis()->SetTitleSize(0.07);
	G_bdtout_data[bdtbin][j][loose]->GetXaxis()->SetNdivisions(5);
	//G_bdtout_data[bdtbin][j][loose]->GetXaxis()->SetTitle("m_{H}");
	G_bdtout_data[bdtbin][j][loose]->Fit("pol1","Q");
	gradient_data[bdtbin][j][loose] = G_bdtout_data[bdtbin][j][loose]->GetFunction("pol1")->GetParameter(1);
	gradient_data_err[bdtbin][j][loose] = G_bdtout_data[bdtbin][j][loose]->GetFunction("pol1")->GetParError(1);

	canvas[j]->cd(bdtbin+1+(nbins[j]*(loose+3)));
	G_bdtout_mc[bdtbin][j][loose] = new TGraphErrors(15,mass,bdtout_mc,mass_err,bdtout_mc_err);
	G_bdtout_mc[bdtbin][j][loose]->SetMarkerSize(.8);
	G_bdtout_mc[bdtbin][j][loose]->SetMarkerColor(4);
	G_bdtout_mc[bdtbin][j][loose]->SetMarkerStyle(20+loose);
	G_bdtout_mc[bdtbin][j][loose]->Draw("AP");
	G_bdtout_mc[bdtbin][j][loose]->GetXaxis()->SetLabelSize(0.07);
	G_bdtout_mc[bdtbin][j][loose]->GetYaxis()->SetLabelSize(0.07);
	G_bdtout_mc[bdtbin][j][loose]->GetXaxis()->SetTitleSize(0.07);
	G_bdtout_mc[bdtbin][j][loose]->GetXaxis()->SetNdivisions(5);
	//G_bdtout_mc[bdtbin][j][loose]->GetXaxis()->SetTitle("m_{H}");
	G_bdtout_mc[bdtbin][j][loose]->Fit("pol1","Q");
	gradient_mc[bdtbin][j][loose] = G_bdtout_mc[bdtbin][j][loose]->GetFunction("pol1")->GetParameter(1);
	gradient_mc_err[bdtbin][j][loose] = G_bdtout_mc[bdtbin][j][loose]->GetFunction("pol1")->GetParError(1);
      }
    }

    Double_t nentries_bin;

    mG_gradient_data[j] = new TMultiGraph();
    mG_gradient_mc[j] = new TMultiGraph();
    mG_gradient[j] = new TMultiGraph();
    mG_gradient_loose[j] = new TMultiGraph();
    mG_gradient_veryloose[j] = new TMultiGraph();
    mG_bias[j] = new TMultiGraph();
    mG_bias_loose[j] = new TMultiGraph();
    mG_bias_veryloose[j] = new TMultiGraph();
    for (int loose=0; loose<3; loose++) {
      for(int bdtbin=0; bdtbin<nbins[j]; bdtbin++) {
	nentries_bin = hist_data_ref[j][loose]->GetBinContent(bdtbin+1);
	gradient_data_temp[bdtbin]=gradient_data[bdtbin][j][loose]/nentries_bin;
	gradient_mc_temp[bdtbin]=gradient_mc[bdtbin][j][loose]/nentries_bin;
	gradient_data_err_temp[bdtbin]=gradient_data_err[bdtbin][j][loose]/nentries_bin;
	gradient_mc_err_temp[bdtbin]=gradient_mc_err[bdtbin][j][loose]/nentries_bin;
	bias_data[bdtbin]=gradient_data_temp[bdtbin]*bias_sf[loose];
	bias_mc[bdtbin]=gradient_mc_temp[bdtbin]*bias_sf[loose];
	bias_data_err[bdtbin]=gradient_data_err_temp[bdtbin]*fabs(bias_sf[loose]);
	bias_mc_err[bdtbin]=gradient_mc_err_temp[bdtbin]*fabs(bias_sf[loose]);
      }
      for(int bdtbin=0; bdtbin<nbins[j]; bdtbin++) bin_shifted[bdtbin]=bin[bdtbin]+binshift[1];
      G_gradient_data[j][loose] = new TGraphErrors(nbins[j],bin_shifted,gradient_data_temp,bin_err,gradient_data_err_temp);
      G_gradient_data[j][loose]->SetMarkerSize(1.);
      G_gradient_data[j][loose]->SetMarkerColor(1);
      G_gradient_data[j][loose]->SetMarkerStyle(20+loose);
      mG_gradient_data[j]->Add(G_gradient_data[j][loose]);
      G_bias_data[j][loose] = new TGraphErrors(nbins[j],bin_shifted,bias_data,bin_err,bias_data_err);
      G_bias_data[j][loose]->SetMarkerSize(1.);
      G_bias_data[j][loose]->SetMarkerColor(1);
      G_bias_data[j][loose]->SetMarkerStyle(20+loose);
      for(int bdtbin=0; bdtbin<nbins[j]; bdtbin++) bin_shifted[bdtbin]=bin[bdtbin]+binshift[2];
      G_gradient_mc[j][loose] = new TGraphErrors(nbins[j],bin_shifted,gradient_mc_temp,bin_err,gradient_mc_err_temp);
      G_gradient_mc[j][loose]->SetMarkerSize(1.);
      G_gradient_mc[j][loose]->SetMarkerColor(4);
      G_gradient_mc[j][loose]->SetMarkerStyle(20+loose);
      mG_gradient_mc[j]->Add(G_gradient_mc[j][loose]);
      G_bias_mc[j][loose] = new TGraphErrors(nbins[j],bin_shifted,bias_mc,bin_err,bias_mc_err);
      G_bias_mc[j][loose]->SetMarkerSize(1.);
      G_bias_mc[j][loose]->SetMarkerColor(4);
      G_bias_mc[j][loose]->SetMarkerStyle(20+loose);
    }
    mG_gradient[j]->Add(G_gradient_data[j][0]);
    mG_gradient_loose[j]->Add(G_gradient_data[j][1]);
    mG_gradient_veryloose[j]->Add(G_gradient_data[j][2]);
    mG_gradient[j]->Add(G_gradient_mc[j][0]);
    mG_gradient_loose[j]->Add(G_gradient_mc[j][1]);
    mG_gradient_veryloose[j]->Add(G_gradient_mc[j][2]);
    mG_bias[j]->Add(G_bias_data[j][0]);
    mG_bias_loose[j]->Add(G_bias_data[j][1]);
    mG_bias_veryloose[j]->Add(G_bias_data[j][2]);
    mG_bias[j]->Add(G_bias_mc[j][0]);
    mG_bias_loose[j]->Add(G_bias_mc[j][1]);
    mG_bias_veryloose[j]->Add(G_bias_mc[j][2]);
  }

  TCanvas *canvas_gradient = new TCanvas("canvas_gradient","canvas_gradient",1400,850);
  canvas_gradient->SetFillColor(0);
  canvas_gradient->Divide(2,3);

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.07);

  TLegend *leg = new TLegend(.15,.15,.35,.35);
  leg->SetBorderSize(0);
  leg->SetTextSize(.05);
  leg->AddEntry(G_gradient_data[0][0],"Data","LP");
  leg->AddEntry(G_gradient_mc[0][1],"MC","LP");

  for(int j=0; j<2; j++) {

    canvas_gradient->cd(1+j);
    gPad->SetGrid();
    mG_gradient[j]->Draw("AP");
    if (zoom) mG_gradient[j]->SetMaximum(0.25);
    if (zoom) mG_gradient[j]->SetMinimum(-0.25);
    mG_gradient[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_gradient[j]->GetXaxis()->SetLabelSize(0.05);
    mG_gradient[j]->GetYaxis()->SetLabelSize(0.05);
    mG_gradient[j]->GetXaxis()->SetTitleSize(0.05);
    mG_gradient[j]->GetYaxis()->SetTitleSize(0.05);
    mG_gradient[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_gradient[j]->GetYaxis()->SetTitle("Slope");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Gradient");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive");
    }
    txt->Draw();
    line[j]->Draw();    
    leg->Draw();    

    canvas_gradient->cd(3+j);
    gPad->SetGrid();
    mG_gradient_loose[j]->Draw("AP");
    if (zoom) mG_gradient_loose[j]->SetMaximum(0.25);
    if (zoom) mG_gradient_loose[j]->SetMinimum(-0.25);
    mG_gradient_loose[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_gradient_loose[j]->GetXaxis()->SetLabelSize(0.05);
    mG_gradient_loose[j]->GetYaxis()->SetLabelSize(0.05);
    mG_gradient_loose[j]->GetXaxis()->SetTitleSize(0.05);
    mG_gradient_loose[j]->GetYaxis()->SetTitleSize(0.05);
    mG_gradient_loose[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_gradient_loose[j]->GetYaxis()->SetTitle("Slope");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Gradient, loose");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive, loose");
    }
    txt->Draw();
    line[j]->Draw(); 
    leg->Draw();    

    canvas_gradient->cd(5+j);
    gPad->SetGrid();
    mG_gradient_veryloose[j]->Draw("AP");
    if (zoom) mG_gradient_veryloose[j]->SetMaximum(0.25);
    if (zoom) mG_gradient_veryloose[j]->SetMinimum(-0.25);
    mG_gradient_veryloose[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_gradient_veryloose[j]->GetXaxis()->SetLabelSize(0.05);
    mG_gradient_veryloose[j]->GetYaxis()->SetLabelSize(0.05);
    mG_gradient_veryloose[j]->GetXaxis()->SetTitleSize(0.05);
    mG_gradient_veryloose[j]->GetYaxis()->SetTitleSize(0.05);
    mG_gradient_veryloose[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_gradient_veryloose[j]->GetYaxis()->SetTitle("Slope");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Gradient, very loose");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive, very loose");
    }
    txt->Draw();
    line[j]->Draw(); 
    leg->Draw();    

  }

  TCanvas *canvas_bias = new TCanvas("canvas_bias","canvas_bias",1400,850);
  canvas_bias->SetFillColor(0);
  canvas_bias->Divide(2,3);

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.07);

  TLegend *leg = new TLegend(.15,.15,.35,.35);
  leg->SetBorderSize(0);
  leg->SetTextSize(.05);
  leg->AddEntry(G_bias_data[0][0],"Data","LP");
  leg->AddEntry(G_bias_mc[0][1],"MC","LP");

  for(int j=0; j<2; j++) {

    canvas_bias->cd(1+j);
    gPad->SetGrid();
    mG_bias[j]->Draw("AP");
    if (zoom) mG_bias[j]->SetMaximum(0.25);
    if (zoom) mG_bias[j]->SetMinimum(-0.25);
    mG_bias[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_bias[j]->GetXaxis()->SetLabelSize(0.05);
    mG_bias[j]->GetYaxis()->SetLabelSize(0.05);
    mG_bias[j]->GetXaxis()->SetTitleSize(0.05);
    mG_bias[j]->GetYaxis()->SetTitleSize(0.05);
    mG_bias[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_bias[j]->GetYaxis()->SetTitle("Bias");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Bias");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive");
    }
    txt->Draw();
    line[j]->Draw();    
    leg->Draw();    

    canvas_bias->cd(3+j);
    gPad->SetGrid();
    mG_bias_loose[j]->Draw("AP");
    if (zoom) mG_bias_loose[j]->SetMaximum(0.25);
    if (zoom) mG_bias_loose[j]->SetMinimum(-0.25);
    mG_bias_loose[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_bias_loose[j]->GetXaxis()->SetLabelSize(0.05);
    mG_bias_loose[j]->GetYaxis()->SetLabelSize(0.05);
    mG_bias_loose[j]->GetXaxis()->SetTitleSize(0.05);
    mG_bias_loose[j]->GetYaxis()->SetTitleSize(0.05);
    mG_bias_loose[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_bias_loose[j]->GetYaxis()->SetTitle("Bias");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Bias, loose");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive, loose");
    }
    txt->Draw();
    line[j]->Draw(); 
    leg->Draw();    

    canvas_bias->cd(5+j);
    gPad->SetGrid();
    mG_bias_veryloose[j]->Draw("AP");
    if (zoom) mG_bias_veryloose[j]->SetMaximum(0.25);
    if (zoom) mG_bias_veryloose[j]->SetMinimum(-0.25);
    mG_bias_veryloose[j]->GetXaxis()->SetRangeUser(0.,float(nbins[j]));
    mG_bias_veryloose[j]->GetXaxis()->SetLabelSize(0.05);
    mG_bias_veryloose[j]->GetYaxis()->SetLabelSize(0.05);
    mG_bias_veryloose[j]->GetXaxis()->SetTitleSize(0.05);
    mG_bias_veryloose[j]->GetYaxis()->SetTitleSize(0.05);
    mG_bias_veryloose[j]->GetXaxis()->SetTitle("BDT output bin number");
    mG_bias_veryloose[j]->GetYaxis()->SetTitle("Bias");
    if (j==0) {
      txt->DrawLatex(0.15,0.4,"Gradient, very loose");
    } else {
      txt->DrawLatex(0.15,0.4,"Adaptive, very loose");
    }
    txt->Draw();
    line[j]->Draw(); 
    leg->Draw();    

  }

  canvas[0]->SaveAs(outdir+"slope_grad_"+mass_str+".gif");
  canvas[1]->SaveAs(outdir+"slope_ada_"+mass_str+".gif");
  if (zoom) {
    canvas_gradient->SaveAs(outdir+"slope_summary_zoomed"+mass_str+".gif");
  } else {
    canvas_gradient->SaveAs(outdir+"slope_summary_"+mass_str+".gif");
  }
  canvas_bias->SaveAs(outdir+"bias_"+mass_str+".gif");

  for(int bdtbin=0; bdtbin<nbins[0]; bdtbin++) {
    nentries_bin = hist_data_ref[0][0]->GetBinContent(bdtbin+1);
    cout << "  slope[" << bdtbin << "] = " << gradient_data[bdtbin][0][0]/nentries_bin << ";" << endl;
  }

}
