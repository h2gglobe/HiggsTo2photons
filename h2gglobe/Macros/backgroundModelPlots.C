#include "TMath.h"
#include <iomanip>

// standard includes
#include <cmath>
#include <map>
#include <set>
#include <vector>

void backgroundModelPlots(int mass_in=120, bool www=false, TString outdirname="BDTplots_all", int sbwidth=2, int loose=0, bool sob=1) {

  bool rebin=true;

  bool madgraph=true;
  bool fakes=true;
  bool data=true;

  bool logy=true;
  float sidebandWidth, nperbin;
  int nSB;
  TString sb_str;
  TString loose_str;
  TString sob_str="";
  if (sob) sob_str="_sob";
  if (sbwidth==2) {
    sidebandWidth=0.02;
    nSB=3;
    sb_str="_2pc";
    if (loose==0) {
      loose_str="";
      nperbin=50.;
    } else if (loose==1) {
      loose_str="_loose";
      nperbin=100.;
    } else if (loose==2) {
      loose_str="_veryloose";
      nperbin=150.;
    } else {
      cout << "loose must be 0, 1 or 2" << endl;
    }
  } else if (sbwidth==7) {
    sidebandWidth=0.07;
    nSB=1;
    sb_str="";
    if (loose>0 && sob) logy=true;
    if (loose==0) {
      loose_str="";
      nperbin=150.;
    } else if (loose==1) {
      loose_str="_loose";
      nperbin=300.;
    } else if (loose==2) {
      loose_str="_veryloose";
      nperbin=450.;
    } else {
      cout << "loose must be 0, 1 or 2" << endl;
    }
  } else {
    cout << "sbwidth must be 2 or 7" << endl;
    break;
  }

  float signalRegionWidth=0.07;
  const float dataMC_sf = 1.;

  bool equalBinWidths=1;
  bool rebinBdtOut=0;

  TString mass_str;
  mass_str+=mass_in;
  TString mass_str2= mass_str;
  if (mass_in==150) mass_str2="145";

  TString cat_str;
  if (mass_in==115) cat_str="0";
  if (mass_in==120) cat_str="1";
  if (mass_in==125) cat_str="2";
  if (mass_in==130) cat_str="3";
  if (mass_in==135) cat_str="4";
  if (mass_in==140) cat_str="5";
  if (mass_in==150) cat_str="6";

  TString outdir;
  if (www) {
    outdir = "/afs/cern.ch/user/f/futyand/www/hgg/"+outdirname+"/"+mass_str+"/gifs/";
  } else {
    outdir = outdirname+"/";
  }

  int sb_low=2;
  int sb_high=4;
  if (nSB>1) {
    sb_low=1;
    sb_high=5;
  }
  if (nSB>2) {
    sb_low=0;
    sb_high=6;
  }

  int nbins_nominal_grad;
  int nbins_nominal_ada;
  Double_t xbins_nominal_grad[100];
  Double_t xbins_nominal_ada[100];
  TString nominal_str="";

  if (rebin) {

    nominal_str="_nominal";

    TH1* hist_nominal_grad;
    TH1* hist_nominal_ada;

    TFile *f_bdtout_nominal = TFile::Open("CMS-HGG_4721pb_9Dec_binning.root");
    //if (sbwidth==2) {
    //  TFile *f_bdtout_nominal = TFile::Open("/tmp/futyand/CMS-HGG_4700pb_02-12-11.root");
    //} else if (sbwidth==7) {
    //  TFile *f_bdtout_nominal = TFile::Open("/tmp/futyand/CMS-HGG_4700pb_02-12-11_7pSB.root");
    //}
    hist_nominal_grad = (TH1*)(f_bdtout_nominal->Get("th1f_data_grad_"+mass_str+".0_cat0"))->Clone();
    hist_nominal_ada = (TH1*)(f_bdtout_nominal->Get("th1f_data_ada_"+mass_str+".0_cat0"))->Clone();

    nbins_nominal_grad = hist_nominal_grad->GetNbinsX();
    for (int ibin=0; ibin<nbins_nominal_grad+1; ibin++) {
      xbins_nominal_grad[ibin] = hist_nominal_grad->GetBinLowEdge(ibin+1);
    }
    nbins_nominal_ada = hist_nominal_ada->GetNbinsX();
    for (int ibin=0; ibin<nbins_nominal_ada+1; ibin++) {
      xbins_nominal_ada[ibin] = hist_nominal_ada->GetBinLowEdge(ibin+1);
    }
    f_bdtout_nominal->Close();

  }

  TFile *f_bdtout = TFile::Open("CMS-HGG_4721pb_9Dec"+loose_str+".root");
  TFile *f_bdtin = TFile::Open("histograms_CMS-HGG_4721pb_9Dec"+loose_str+".root");

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  TH1* hist_sig[24];
  TH1* hist_data[8][24];
  TH1* hist_born[8][24];
  TH1* hist_box[8][24];
  TH1* hist_gjet_pp[8][24];
  TH1* hist_gjet_pf[8][24];
  TH1* hist_qcd_pp[8][24];
  TH1* hist_qcd_pf[8][24];
  TH1* hist_qcd_ff[8][24];
  TH1* hist_dy[8][24];
  TH1* hist_bkg[8][24];
  TH1* hist_bkgModel[24];
  TH1* hist_balg[24];
  TH1* hist_bkg_scaled[8][24];
  TH1* hist_data_scaled[8][24];
  TH1* hist_bkg_sig[24];
  TH1* hist_data_sig[24];
  TH1* hist_sig_reweight[24];
  TH1* hist_bkg_reweight[24];
  TH1* hist_data_reweight[24];
  TH1* hist_data_reweight_rebin[24];
  TH1* hist_bkg_ratio[24];
  TH1* hist_bkg_ratio_low[24];
  TH1* hist_bkg_ratio_high[24];
  TH1* hist_data_ratio[24];
  THStack* hist_bkg_stack[24];
  THStack* hist_bkg_stack_sig[24];

  TH1 *hist_mass_sig;
  TH1 *hist_mass_data;
  TH1 *hist_mass_born;
  TH1 *hist_mass_box;
  TH1 *hist_mass_gjet_pp;
  TH1 *hist_mass_gjet_pf;
  TH1 *hist_mass_qcd_pp;
  TH1 *hist_mass_qcd_pf;
  TH1 *hist_mass_qcd_ff;
  TH1 *hist_mass_dy;
  TH1 *hist_mass_bkg;
  THStack* hist_mass_bkg_stack;

  TH1* hist_sig_fine[24];
  TH1* hist_data_fine[8][24];
  TH1* hist_bkg_fine[8][24];
  TH1* hist_bkgModel_fine[24];
  TH1* hist_balg_fine[24];

  f_bdtout->cd();

  if (!equalBinWidths) {

    hist_sig[22] = (TH1*)(f_bdtout->Get("th1f_sig_grad_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_data[0][22] = (TH1*)(f_bdtout->Get("th1f_bkg_3low_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_data[1][22] = (TH1*)(f_bdtout->Get("th1f_bkg_2low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data[2][22] = (TH1*)(f_bdtout->Get("th1f_bkg_1low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data[3][22] = (TH1*)(f_bdtout->Get("th1f_data_grad_"+mass_str+".0_cat0"))->Clone();
    hist_data[4][22] = (TH1*)(f_bdtout->Get("th1f_bkg_1high_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_data[5][22] = (TH1*)(f_bdtout->Get("th1f_bkg_2high_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_data[6][22] = (TH1*)(f_bdtout->Get("th1f_bkg_3high_grad_"+mass_str+".0_cat0"))->Clone();
    hist_bkgModel[22] = (TH1*)(f_bdtout->Get("th1f_bkg_grad_"+mass_str+".0_cat0"))->Clone();
    hist_balg[22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_balg_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_bkg[0][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_3low_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_bkg[1][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_2low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[2][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_1low_grad_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[3][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_grad_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[4][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_1high_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_bkg[5][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_2high_grad_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_bkg[6][22] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_3high_grad_"+mass_str+".0_cat0"))->Clone();

    hist_sig[23] = (TH1*)(f_bdtout->Get("th1f_sig_ada_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_data[0][23] = (TH1*)(f_bdtout->Get("th1f_bkg_3low_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_data[1][23] = (TH1*)(f_bdtout->Get("th1f_bkg_2low_ada_"+mass_str+".0_cat0"))->Clone();
    hist_data[2][23] = (TH1*)(f_bdtout->Get("th1f_bkg_1low_ada_"+mass_str+".0_cat0"))->Clone();
    hist_data[3][23] = (TH1*)(f_bdtout->Get("th1f_data_ada_"+mass_str+".0_cat0"))->Clone();
    hist_data[4][23] = (TH1*)(f_bdtout->Get("th1f_bkg_1high_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_data[5][23] = (TH1*)(f_bdtout->Get("th1f_bkg_2high_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_data[6][23] = (TH1*)(f_bdtout->Get("th1f_bkg_3high_ada_"+mass_str+".0_cat0"))->Clone();
    hist_bkgModel[23] = (TH1*)(f_bdtout->Get("th1f_bkg_ada_"+mass_str+".0_cat0"))->Clone();
    hist_balg[23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_balg_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_bkg[0][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_3low_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_bkg[1][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_2low_ada_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[2][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_1low_ada_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[3][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_ada_"+mass_str+".0_cat0"))->Clone();
    hist_bkg[4][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_1high_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>1) hist_bkg[5][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_2high_ada_"+mass_str+".0_cat0"))->Clone();
    if (nSB>2) hist_bkg[6][23] = (TH1*)(f_bdtout->Get("th1f_bkg_mc_3high_ada_"+mass_str+".0_cat0"))->Clone();

  } else {

    if (rebin) {

      hist_sig_fine[22] = (TH1*)(f_bdtout->Get("th1f_sig_BDT_grad_ggh_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_data_fine[0][22] = (TH1*)(f_bdtout->Get("th1f_data_3low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_data_fine[1][22] = (TH1*)(f_bdtout->Get("th1f_data_2low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[2][22] = (TH1*)(f_bdtout->Get("th1f_data_1low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[3][22] = (TH1*)(f_bdtout->Get("th1f_data_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[4][22] = (TH1*)(f_bdtout->Get("th1f_data_1high_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_data_fine[5][22] = (TH1*)(f_bdtout->Get("th1f_data_2high_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_data_fine[6][22] = (TH1*)(f_bdtout->Get("th1f_data_3high_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_bkgModel_fine[22] = (TH1*)(f_bdtout->Get("th1f_data_BDT_sideband_grad_"+mass_str+".0_cat0"))->Clone();
      hist_balg_fine[22] = (TH1*)(f_bdtout->Get("th1f_bkg_BDT_grad_all_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_bkg_fine[0][22] = (TH1*)(f_bdtout->Get("th1f_bkg_3low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_bkg_fine[1][22] = (TH1*)(f_bdtout->Get("th1f_bkg_2low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[2][22] = (TH1*)(f_bdtout->Get("th1f_bkg_1low_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[3][22] = (TH1*)(f_bdtout->Get("th1f_bkg_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[4][22] = (TH1*)(f_bdtout->Get("th1f_bkg_1high_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_bkg_fine[5][22] = (TH1*)(f_bdtout->Get("th1f_bkg_2high_BDT_grad_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_bkg_fine[6][22] = (TH1*)(f_bdtout->Get("th1f_bkg_3high_BDT_grad_"+mass_str+".0_cat0"))->Clone();

      hist_sig_fine[23] = (TH1*)(f_bdtout->Get("th1f_sig_BDT_ada_ggh_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_data_fine[0][23] = (TH1*)(f_bdtout->Get("th1f_data_3low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_data_fine[1][23] = (TH1*)(f_bdtout->Get("th1f_data_2low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[2][23] = (TH1*)(f_bdtout->Get("th1f_data_1low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[3][23] = (TH1*)(f_bdtout->Get("th1f_data_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_data_fine[4][23] = (TH1*)(f_bdtout->Get("th1f_data_1high_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_data_fine[5][23] = (TH1*)(f_bdtout->Get("th1f_data_2high_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_data_fine[6][23] = (TH1*)(f_bdtout->Get("th1f_data_3high_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_bkgModel_fine[23] = (TH1*)(f_bdtout->Get("th1f_data_BDT_sideband_ada_"+mass_str+".0_cat0"))->Clone();
      hist_balg_fine[22] = (TH1*)(f_bdtout->Get("th1f_bkg_BDT_ada_all_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_bkg_fine[0][23] = (TH1*)(f_bdtout->Get("th1f_bkg_3low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_bkg_fine[1][23] = (TH1*)(f_bdtout->Get("th1f_bkg_2low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[2][23] = (TH1*)(f_bdtout->Get("th1f_bkg_1low_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[3][23] = (TH1*)(f_bdtout->Get("th1f_bkg_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      hist_bkg_fine[4][23] = (TH1*)(f_bdtout->Get("th1f_bkg_1high_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>1) hist_bkg_fine[5][23] = (TH1*)(f_bdtout->Get("th1f_bkg_2high_BDT_ada_"+mass_str+".0_cat0"))->Clone();
      if (nSB>2) hist_bkg_fine[6][23] = (TH1*)(f_bdtout->Get("th1f_bkg_3high_BDT_ada_"+mass_str+".0_cat0"))->Clone();

      hist_sig_fine[22]->Rebin(nbins_nominal_grad,"th1f_nominal_sig_grad_ggh_"+mass_str+".0_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>2) hist_data_fine[0][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_3low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>1) hist_data_fine[1][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_2low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_data_fine[2][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_1low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_data_fine[3][22]->Rebin(nbins_nominal_grad,"th1f_nominal_data_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_data_fine[4][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_1high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>1) hist_data_fine[5][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_2high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>2) hist_data_fine[6][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_3high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_bkgModel_fine[22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_balg_fine[22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_balg_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>2) hist_bkg_fine[0][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_3low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>1) hist_bkg_fine[1][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_2low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_bkg_fine[2][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_1low_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_bkg_fine[3][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      hist_bkg_fine[4][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_1high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>1) hist_bkg_fine[5][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_2high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);
      if (nSB>2) hist_bkg_fine[6][22]->Rebin(nbins_nominal_grad,"th1f_nominal_bkg_mc_3high_grad_"+mass_str+".0_cat0",xbins_nominal_grad);

      hist_sig_fine[23]->Rebin(nbins_nominal_ada,"th1f_nominal_sig_ada_ggh_"+mass_str+".0_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>2) hist_data_fine[0][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_3low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>1) hist_data_fine[1][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_2low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_data_fine[2][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_1low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_data_fine[3][23]->Rebin(nbins_nominal_ada,"th1f_nominal_data_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_data_fine[4][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_1high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>1) hist_data_fine[5][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_2high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>2) hist_data_fine[6][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_3high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_bkgModel_fine[23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_balg_fine[22]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_balg_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>2) hist_bkg_fine[0][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_3low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>1) hist_bkg_fine[1][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_2low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_bkg_fine[2][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_1low_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_bkg_fine[3][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      hist_bkg_fine[4][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_1high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>1) hist_bkg_fine[5][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_2high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);
      if (nSB>2) hist_bkg_fine[6][23]->Rebin(nbins_nominal_ada,"th1f_nominal_bkg_mc_3high_ada_"+mass_str+".0_cat0",xbins_nominal_ada);

    }

    int nbins = ((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_grad_"+mass_str+".0_cat0"))->GetNbinsX();
    int binshift=0;
    if (rebinBdtOut && nbins%2!=0) {
      nbins-=1;
      binshift=1;
    }
    hist_sig[22] = new TH1F("hist_sig_grad_ggh","hist_sig_grad_ggh",nbins,0,float(nbins));
    if (nSB>2) hist_data[0][22] = new TH1F("hist_data_3low_grad","hist_data_3low_grad",nbins,0,float(nbins));
    if (nSB>1) hist_data[1][22] = new TH1F("hist_data_2low_grad","hist_data_2low_grad",nbins,0,float(nbins));
    hist_data[2][22] = new TH1F("hist_data_1low_grad","hist_data_1low_grad",nbins,0,float(nbins));
    hist_data[3][22] = new TH1F("hist_data_sig_grad_ggh","hist_data_sig_grad_ggh",nbins,0,float(nbins));
    hist_data[4][22] = new TH1F("hist_data_1high_grad","hist_data_1high_grad",nbins,0,float(nbins));
    if (nSB>1) hist_data[5][22] = new TH1F("hist_data_2high_grad","hist_data_2high_grad",nbins,0,float(nbins));
    if (nSB>2) hist_data[6][22] = new TH1F("hist_data_3high_grad","hist_data_3high_grad",nbins,0,float(nbins));
    hist_bkgModel[22] = new TH1F("hist_bkgModel_grad","hist_bkgModel_grad",nbins,0,float(nbins));
    hist_balg[22] = new TH1F("hist_balg_grad","hist_balg_grad",nbins,0,float(nbins));
    if (nSB>2) hist_bkg[0][22] = new TH1F("hist_bkg_3low_grad","hist_bkg_3low_grad",nbins,0,float(nbins));
    if (nSB>1) hist_bkg[1][22] = new TH1F("hist_bkg_2low_grad","hist_bkg_2low_grad",nbins,0,float(nbins));
    hist_bkg[2][22] = new TH1F("hist_bkg_1low_grad","hist_bkg_1low_grad",nbins,0,float(nbins));
    hist_bkg[3][22] = new TH1F("hist_bkg_sig_grad_ggh","hist_bkg_sig_grad_ggh",nbins,0,float(nbins));
    hist_bkg[4][22] = new TH1F("hist_bkg_1high_grad","hist_bkg_1high_grad",nbins,0,float(nbins));
    if (nSB>1) hist_bkg[5][22] = new TH1F("hist_bkg_2high_grad","hist_bkg_2high_grad",nbins,0,float(nbins));
    if (nSB>2) hist_bkg[6][22] = new TH1F("hist_bkg_3high_grad","hist_bkg_3high_grad",nbins,0,float(nbins));
    for (int ibin=0; ibin<nbins+1; ibin++) {
      hist_sig[22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_sig_grad_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_sig[22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_sig_grad_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_data[0][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_data[0][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_data[1][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_data[1][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[2][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[2][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[3][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[3][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[4][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[4][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_data[5][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_data[5][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_data[6][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_data[6][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkgModel[22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkgModel[22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_balg[22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_balg_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_balg[22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_balg_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[0][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[0][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_bkg[1][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_bkg[1][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[2][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1low_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[2][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1low_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[3][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[3][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[4][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[4][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_bkg[5][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_bkg[5][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_bkg[6][22]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3high_grad_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[6][22]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3high_grad_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
    }

    nbins = ((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_ada_"+mass_str+".0_cat0"))->GetNbinsX();
    binshift=0;
    if (rebinBdtOut && nbins%2!=0) {
      nbins-=1;
      binshift=1;
    }
    hist_sig[23] = new TH1F("hist_sig_ada_ggh","hist_sig_ada_ggh",nbins,0,float(nbins));
    if (nSB>2) hist_data[0][23] = new TH1F("hist_data_3low_ada","hist_data_3low_ada",nbins,0,float(nbins));
    if (nSB>1) hist_data[1][23] = new TH1F("hist_data_2low_ada","hist_data_2low_ada",nbins,0,float(nbins));
    hist_data[2][23] = new TH1F("hist_data_1low_ada","hist_data_1low_ada",nbins,0,float(nbins));
    hist_data[3][23] = new TH1F("hist_data_sig_ada_ggh","hist_data_sig_ada_ggh",nbins,0,float(nbins));
    hist_data[4][23] = new TH1F("hist_data_1high_ada","hist_data_1high_ada",nbins,0,float(nbins));
    if (nSB>1) hist_data[5][23] = new TH1F("hist_data_2high_ada","hist_data_2high_ada",nbins,0,float(nbins));
    if (nSB>2) hist_data[6][23] = new TH1F("hist_data_3high_ada","hist_data_3high_ada",nbins,0,float(nbins));
    hist_bkgModel[23] = new TH1F("hist_bkgModel_ada","hist_bkgModel_ada",nbins,0,float(nbins));
    hist_balg[23] = new TH1F("hist_balg_ada","hist_balg_ada",nbins,0,float(nbins));
    if (nSB>2) hist_bkg[0][23] = new TH1F("hist_bkg_3low_ada","hist_bkg_3low_ada",nbins,0,float(nbins));
    if (nSB>1) hist_bkg[1][23] = new TH1F("hist_bkg_2low_ada","hist_bkg_2low_ada",nbins,0,float(nbins));
    hist_bkg[2][23] = new TH1F("hist_bkg_1low_ada","hist_bkg_1low_ada",nbins,0,float(nbins));
    hist_bkg[3][23] = new TH1F("hist_bkg_sig_ada_ggh","hist_bkg_sig_ada_ggh",nbins,0,float(nbins));
    hist_bkg[4][23] = new TH1F("hist_bkg_1high_ada","hist_bkg_1high_ada",nbins,0,float(nbins));
    if (nSB>1) hist_bkg[5][23] = new TH1F("hist_bkg_2high_ada","hist_bkg_2high_ada",nbins,0,float(nbins));
    if (nSB>2) hist_bkg[6][23] = new TH1F("hist_bkg_3high_ada","hist_bkg_3high_ada",nbins,0,float(nbins));
    for (int ibin=0; ibin<nbins+1; ibin++) {
      hist_sig[23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_sig_ada_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_sig[23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_sig_ada_ggh_"+mass_str+".0_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_data[0][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_data[0][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_data[1][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_data[1][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[2][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[2][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[3][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[3][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_data_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_data[4][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_data[4][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_1high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_data[5][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_data[5][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_2high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_data[6][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_data[6][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_3high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkgModel[23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkgModel[23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_balg[23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_balg_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_balg[23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_balg_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[0][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[0][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_bkg[1][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_bkg[1][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[2][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1low_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[2][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1low_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[3][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[3][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      hist_bkg[4][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      hist_bkg[4][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_1high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>1) hist_bkg[5][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>1) hist_bkg[5][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_2high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
      if (nSB>2) hist_bkg[6][23]->SetBinContent(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3high_ada_"+mass_str+".0_cat0"))->GetBinContent(ibin+binshift));
      if (nSB>2) hist_bkg[6][23]->SetBinError(ibin,((TH1*)f_bdtout->Get("th1f"+nominal_str+"_bkg_mc_3high_ada_"+mass_str+".0_cat0"))->GetBinError(ibin+binshift));
    }

  }

  for (int j=22; j<24; j++) {
    if (nSB<3) {
      hist_data[0][j]=(TH1*)hist_data[3][j]->Clone();
      hist_data[6][j]=(TH1*)hist_data[3][j]->Clone();
      hist_bkg[0][j]=(TH1*)hist_bkg[3][j]->Clone();
      hist_bkg[6][j]=(TH1*)hist_bkg[3][j]->Clone();
      hist_data[0][j]->Reset();
      hist_data[6][j]->Reset();
      hist_bkg[0][j]->Reset();
      hist_bkg[6][j]->Reset();
    }
    if (nSB<2) {
      hist_data[1][j]=(TH1*)hist_data[3][j]->Clone();
      hist_data[5][j]=(TH1*)hist_data[3][j]->Clone();
      hist_bkg[1][j]=(TH1*)hist_bkg[3][j]->Clone();
      hist_bkg[5][j]=(TH1*)hist_bkg[3][j]->Clone();
      hist_data[1][j]->Reset();
      hist_data[5][j]->Reset();
      hist_bkg[1][j]->Reset();
      hist_bkg[5][j]->Reset();
    }
    if (!equalBinWidths) {
      hist_sig[j]->GetYaxis()->SetTitle("");
      hist_bkgModel[j]->GetYaxis()->SetTitle("");
      hist_balg[j]->GetYaxis()->SetTitle("");
      for (int i=0; i<7; i++) {
	hist_data[i][j]->GetYaxis()->SetTitle("");
	hist_bkg[i][j]->GetYaxis()->SetTitle("");
      }
    }

    if (rebinBdtOut) {
      hist_sig[j]->Rebin(2);
      for (int i=0; i<7; i++) {
	hist_data[i][j]->Rebin(2);
	hist_bkg[i][j]->Rebin(2);
	hist_bkgModel[j]->Rebin(2);
	hist_balg[j]->Rebin(2);
      }
    }

  }

  f_bdtin->cd();

  hist_mass_sig = (TH1*)(f_bdtin->Get("all_mass_cat0_gluglu_H_gg_"+mass_str2+"_pu2011"));
  if (data) hist_mass_data = (TH1*)all_mass_cat0_Data->Clone();
  if (!madgraph) {
    hist_mass_born = (TH1*)all_mass_cat0_Born25->Clone();
    hist_mass_gjet_pp = (TH1*)all_mass_cat0_GJetPP->Clone();
    hist_mass_qcd_pp = (TH1*)all_mass_cat0_QCDPP->Clone();
  } else {
    hist_mass_born = (TH1*)all_mass_cat0_DiPhotonJets->Clone();
  }
  hist_mass_box = (TH1*)all_mass_cat0_Box25->Clone();
  hist_mass_gjet_pf = (TH1*)all_mass_cat0_GJet->Clone();
  hist_mass_qcd_pf = (TH1*)all_mass_cat0_QCD40PF->Clone();
  hist_mass_qcd_pf->Add(all_mass_cat0_QCD30PF);
  hist_mass_qcd_ff = (TH1*)all_mass_cat0_QCD40FF->Clone();
  hist_mass_qcd_ff->Add(all_mass_cat0_QCD30FF);
  hist_mass_dy = (TH1*)all_mass_cat0_DYJetsToLL->Clone();

  TString var[24] = {"ptOverMH","eta","deltaPhi","cosDeltaPhi","pho1_pt","pho2_pt","pho1_eta","pho2_eta","pho_minr9","maxeta","ptOverMH","pho1_ptOverMH","pho2_ptOverMH","sigmaMOverM","deltaEta","deltaMOverMH","pho1_ptOverMH","pho2_ptOverMH","sigmaMOverMH","sigmaMOverMH_Eonly","deltaMOverSigmaM","sigmaM","bdtOut_grad","bdtOut_ada"};

  for (int ivar=0; ivar<1; ivar++) {
    if (ivar==10 || ivar==16 || ivar==17) continue;
    cout << var[ivar] << endl;
    hist_sig[ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_gluglu_H_gg_"+mass_str2+"_pu2011"))->Clone();

    if (data) {
      hist_data[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_Data"))->Clone();
      hist_data[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_Data"))->Clone();
      hist_data[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_Data"))->Clone();
      hist_data[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_Data"))->Clone();
      hist_data[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_Data"))->Clone();
      hist_data[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_Data"))->Clone();
      hist_data[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_Data"))->Clone();
    }

    if (!madgraph) {

      hist_born[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_Born25"))->Clone();
      hist_born[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_Born25"))->Clone();
      hist_born[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_Born25"))->Clone();
      hist_born[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_Born25"))->Clone();
      hist_born[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_Born25"))->Clone();
      hist_born[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_Born25"))->Clone();
      hist_born[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_Born25"))->Clone();

      hist_gjet_pp[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_GJetPP"))->Clone();
      hist_gjet_pp[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_GJetPP"))->Clone();

      hist_qcd_pp[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_QCDPP"))->Clone();
      hist_qcd_pp[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_QCDPP"))->Clone();

    } else {

      hist_born[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_DiPhotonJets"))->Clone();
      hist_born[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_DiPhotonJets"))->Clone();

    }

    hist_box[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_Box25"))->Clone();
    hist_box[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_Box25"))->Clone();
    hist_box[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_Box25"))->Clone();
    hist_box[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_Box25"))->Clone();
    hist_box[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_Box25"))->Clone();
    hist_box[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_Box25"))->Clone();
    hist_box[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_Box25"))->Clone();

    hist_gjet_pf[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_GJet"))->Clone();
    hist_gjet_pf[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_GJet"))->Clone();

    hist_qcd_pf[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_QCD40PF"))->Clone();
    hist_qcd_pf[0][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[1][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[2][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[3][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[4][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[5][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_QCD30PF")));
    hist_qcd_pf[6][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_QCD30PF")));

    hist_qcd_ff[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_QCD40FF"))->Clone();
    hist_qcd_ff[0][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[1][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[2][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[3][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[4][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[5][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_QCD30FF")));
    hist_qcd_ff[6][ivar]->Add((TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_QCD30FF")));

    hist_dy[0][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow3_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[1][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow2_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[2][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mlow_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[3][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_msig_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[4][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[5][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh2_cat"+cat_str+"_DYJetsToLL"))->Clone();
    hist_dy[6][ivar] = (TH1*)(f_bdtin->Get(var[ivar]+"_mhigh3_cat"+cat_str+"_DYJetsToLL"))->Clone();

  }

  TString title[24] = {"Diphoton p_{T}/M_{H}","Diphoton #eta","#Delta#phi","cos#Delta#phi","lead p_{T} (GeV)","sublead p_{T} (GeV)","lead #eta","sublead #eta","min(R9)","max #eta","Diphoton p_{T}/M_{H}","p_{T1}/M_{H}","p_{T2}/M_{H}","#sigma_{M}/M_{#gamma#gamma}","#Delta#eta","#DeltaM/M_{H}","p_{T1}/M_{H}","p_{T2}/M_{H}","#sigma_{M}/M_{H}","#sigma_{M}/M_{H} (energy only)","#DeltaM/#sigma_{M}","#sigma_{M}","BDT output (gradient boost)","BDT output (adaptive boost)"};
  if (equalBinWidths) {
    title[22]="BDT output bin number (gradient boost)";
    title[23]="BDT output bin number (adaptive boost)";
    var[22]="bdtOutBin_grad";
    var[23]="bdtOutBin_ada";
  }

  TCanvas* canvas[24];

  //TCanvas *c_mgg2 = new TCanvas("c_mgg2","Mgg graph");
  //c_mgg2->SetFillColor(0);

  float mass_hypothesis = float(mass_in);

  /*
  float sidebandShift = 0.;

  double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth) - sidebandShift;
  float mass_hypothesis_low1 = (mass_hypothesis*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,0)) - sidebandShift;
  float mass_hypothesis_low2 = (mass_hypothesis*(1.-sidebandWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,1)) - sidebandShift;
  float mass_hypothesis_low3 = (mass_hypothesis*(1.-sidebandWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,2)) - sidebandShift;

  hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth) + sidebandShift;
  float mass_hypothesis_high1 = (mass_hypothesis*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,0)) + sidebandShift;
  float mass_hypothesis_high2 = (mass_hypothesis*(1.+sidebandWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,1)) + sidebandShift;
  float mass_hypothesis_high3 = (mass_hypothesis*(1.+sidebandWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,2)) + sidebandShift;
  */

  float mass_hypothesis_low1 = mass_hypothesis*(1-signalRegionWidth)/(1+sidebandWidth);
  float mass_hypothesis_low2 = mass_hypothesis_low1*(1-sidebandWidth)/(1+sidebandWidth);
  float mass_hypothesis_low3 = mass_hypothesis_low2*(1-sidebandWidth)/(1+sidebandWidth);
  float mass_hypothesis_high1 = mass_hypothesis*(1+signalRegionWidth)/(1-sidebandWidth);
  float mass_hypothesis_high2 = mass_hypothesis_high1*(1+sidebandWidth)/(1-sidebandWidth);
  float mass_hypothesis_high3 = mass_hypothesis_high2*(1+sidebandWidth)/(1-sidebandWidth);

  //define sidebands
  float sideband_boundaries[8];
  sideband_boundaries[0] = mass_hypothesis_low3*(1-sidebandWidth);
  sideband_boundaries[1] = mass_hypothesis_low2*(1-sidebandWidth);
  sideband_boundaries[2] = mass_hypothesis_low1*(1-sidebandWidth);
  sideband_boundaries[3] = mass_hypothesis*(1-signalRegionWidth);
  sideband_boundaries[4] = mass_hypothesis*(1+signalRegionWidth);
  sideband_boundaries[5] = mass_hypothesis_high1*(1+sidebandWidth);
  sideband_boundaries[6] = mass_hypothesis_high2*(1+sidebandWidth);
  sideband_boundaries[7] = mass_hypothesis_high3*(1+sidebandWidth);

  float signalRegion_boundaries[2];
  signalRegion_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
  signalRegion_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

  for (int i=0; i<8; i++) cout << sideband_boundaries[i] << " ";
  cout << endl;

  float frac_signal = hist_mass_sig->Integral(hist_mass_sig->FindBin(signalRegion_boundaries[0]),hist_mass_sig->FindBin(signalRegion_boundaries[1]))/hist_mass_sig->Integral();
  cout << "fraction of signal in signal region = " << frac_signal << endl;

  float Ndata_sig = hist_bkgModel[22]->Integral();
  cout << Ndata_sig << endl;

  for (int ivar=0; ivar<1; ivar++) {

    if (ivar==10 || ivar==16 || ivar==17) continue;

    if (ivar==0 || ivar==10) {
      hist_sig[ivar]->Rebin(2);
      for (int i=0; i<7; i++) {
	hist_data[i][ivar]->Rebin(2);
	hist_born[i][ivar]->Rebin(2);
	hist_box[i][ivar]->Rebin(2);
	hist_gjet_pf[i][ivar]->Rebin(2);
	hist_qcd_pf[i][ivar]->Rebin(2);
 	hist_qcd_ff[i][ivar]->Rebin(2);
	hist_dy[i][ivar]->Rebin(2);
	if (!madgraph) {
	  hist_gjet_pp[i][ivar]->Rebin(2);
	  hist_qcd_pp[i][ivar]->Rebin(2);
	}
      }
    }

    if (!fakes) {
      for (int i=0; i<7; i++) {
	hist_born[i][ivar]->Scale(1./0.6767);
	hist_box[i][ivar]->Scale(1./0.6767);
	if (!madgraph) {
	  hist_gjet_pp[i][ivar]->Scale(1./0.6767);
	  hist_qcd_pp[i][ivar]->Scale(1./0.6767);
	}
      }
    }

    if (data) hist_data[7][ivar] = (TH1*)hist_data[0][ivar]->Clone();
    hist_born[7][ivar] = (TH1*)hist_born[0][ivar]->Clone();
    hist_box[7][ivar] = (TH1*)hist_box[0][ivar]->Clone();
    hist_gjet_pf[7][ivar] = (TH1*)hist_gjet_pf[0][ivar]->Clone();
    hist_qcd_pf[7][ivar] = (TH1*)hist_qcd_pf[0][ivar]->Clone();
    hist_qcd_ff[7][ivar] = (TH1*)hist_qcd_ff[0][ivar]->Clone();
    hist_dy[7][ivar] = (TH1*)hist_dy[0][ivar]->Clone();
    if (!madgraph) {
      hist_gjet_pp[7][ivar] = (TH1*)hist_gjet_pp[0][ivar]->Clone();
      hist_qcd_pp[7][ivar] = (TH1*)hist_qcd_pp[0][ivar]->Clone();
    }
    for (int i=1; i<7; i++) {
      if (data) hist_data[7][ivar]->Add(hist_data[i][ivar]);
      hist_born[7][ivar]->Add(hist_born[i][ivar]);
      hist_box[7][ivar]->Add(hist_box[i][ivar]);
      hist_gjet_pf[7][ivar]->Add(hist_gjet_pf[i][ivar]);
      hist_qcd_pf[7][ivar]->Add(hist_qcd_pf[i][ivar]);
      hist_qcd_ff[7][ivar]->Add(hist_qcd_ff[i][ivar]);
      hist_dy[7][ivar]->Add(hist_dy[i][ivar]);
      if (!madgraph) {
	hist_gjet_pp[7][ivar]->Add(hist_gjet_pp[i][ivar]);
	hist_qcd_pp[7][ivar]->Add(hist_qcd_pp[i][ivar]);
      }
    }

    /*
    //Apply k-factors
    float bornSF = madgraph ? 1.15 : 1.3;
    for (int i=0; i<8; i++) {
      hist_born[i][ivar]->Scale(bornSF*dataMC_sf);
      hist_box[i][ivar]->Scale(1.3*dataMC_sf);
      hist_gjet_pf[i][ivar]->Scale(1.3*dataMC_sf);
      hist_qcd_pf[i][ivar]->Scale(1.3*dataMC_sf);
      hist_qcd_ff[i][ivar]->Scale(1.*dataMC_sf);
      hist_dy[i][ivar]->Scale(1.15*2321./992.*dataMC_sf);
      if (!madgraph) {
	hist_gjet_pp[i][ivar]->Scale(1.3*dataMC_sf);
	hist_qcd_pp[i][ivar]->Scale(1.3*dataMC_sf);
      }
    }
    */

    for (int i=0; i<8; i++) {
      hist_born[i][ivar]->Scale(dataMC_sf);
      hist_box[i][ivar]->Scale(dataMC_sf);
      hist_gjet_pf[i][ivar]->Scale(dataMC_sf);
      hist_qcd_pf[i][ivar]->Scale(dataMC_sf);
      hist_qcd_ff[i][ivar]->Scale(dataMC_sf);
      hist_dy[i][ivar]->Scale(dataMC_sf);
      if (!madgraph) {
	hist_gjet_pp[i][ivar]->Scale(dataMC_sf);
	hist_qcd_pp[i][ivar]->Scale(dataMC_sf);
      }
    }

    //Apply +ve/-ve folding for eta distributions
    if (ivar==1 || ivar==6 ||ivar==7 || ivar==9) {
      for (int ibin=25; ibin<50; ibin++) {
	hist_sig[ivar]->SetBinContent(ibin+1,hist_sig[ivar]->GetBinContent(ibin+1)+hist_sig[ivar]->GetBinContent(50-ibin));
	for (int i=0; i<8; i++) {
	  hist_data[i][ivar]->SetBinContent(ibin+1,hist_data[i][ivar]->GetBinContent(ibin+1)+hist_data[i][ivar]->GetBinContent(50-ibin));
	  hist_born[i][ivar]->SetBinContent(ibin+1,hist_born[i][ivar]->GetBinContent(ibin+1)+hist_born[i][ivar]->GetBinContent(50-ibin));
	  hist_box[i][ivar]->SetBinContent(ibin+1,hist_box[i][ivar]->GetBinContent(ibin+1)+hist_box[i][ivar]->GetBinContent(50-ibin));
	  hist_gjet_pf[i][ivar]->SetBinContent(ibin+1,hist_gjet_pf[i][ivar]->GetBinContent(ibin+1)+hist_gjet_pf[i][ivar]->GetBinContent(50-ibin));
	  hist_qcd_pf[i][ivar]->SetBinContent(ibin+1,hist_qcd_pf[i][ivar]->GetBinContent(ibin+1)+hist_qcd_pf[i][ivar]->GetBinContent(50-ibin));
	  hist_qcd_ff[i][ivar]->SetBinContent(ibin+1,hist_qcd_ff[i][ivar]->GetBinContent(ibin+1)+hist_qcd_ff[i][ivar]->GetBinContent(50-ibin));
	  hist_dy[i][ivar]->SetBinContent(ibin+1,hist_dy[i][ivar]->GetBinContent(ibin+1)+hist_dy[i][ivar]->GetBinContent(50-ibin));
	  if (!madgraph) {
	    hist_gjet_pp[i][ivar]->SetBinContent(ibin+1,hist_gjet_pp[i][ivar]->GetBinContent(ibin+1)+hist_gjet_pp[i][ivar]->GetBinContent(50-ibin));
	    hist_qcd_pp[i][ivar]->SetBinContent(ibin+1,hist_qcd_pp[i][ivar]->GetBinContent(ibin+1)+hist_qcd_pp[i][ivar]->GetBinContent(50-ibin));
	  }
	}
      }
      for (int ibin=0; ibin<25; ibin++) {
	hist_sig[ivar]->SetBinContent(ibin+1,0.);
	for (int i=0; i<8; i++) {
	  hist_data[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_born[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_box[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_gjet_pf[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_qcd_pf[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_qcd_ff[i][ivar]->SetBinContent(ibin+1,0.);
	  hist_dy[i][ivar]->SetBinContent(ibin+1,0.);
	  if (!madgraph) {
	    hist_gjet_pp[i][ivar]->SetBinContent(ibin+1,0.);
	    hist_qcd_pp[i][ivar]->SetBinContent(ibin+1,0.);
	  }
	}
      }
    }

    for (int i=0; i<8; i++) {
      SetHistogramErrors(hist_data[i][ivar]);
      SetHistogramErrors(hist_born[i][ivar]);
      SetHistogramErrors(hist_box[i][ivar]);
      SetHistogramErrors(hist_gjet_pf[i][ivar]);
      SetHistogramErrors(hist_qcd_pf[i][ivar]);
      SetHistogramErrors(hist_qcd_ff[i][ivar]);
      SetHistogramErrors(hist_dy[i][ivar]);
      if (!madgraph) {
	SetHistogramErrors(hist_gjet_pp[i][ivar]);
	SetHistogramErrors(hist_qcd_pp[i][ivar]);
      }
    }

    hist_sig[ivar]->SetLineColor(4);
    hist_sig[ivar]->SetLineWidth(2.5);
    if (data) hist_data[7][ivar]->SetMarkerStyle(20);
    if (data) hist_data[7][ivar]->SetMarkerSize(.8);
    if (data) hist_data[3][ivar]->SetMarkerStyle(20);
    if (data) hist_data[3][ivar]->SetMarkerSize(.8);
    for (int i=0; i<8; i++) {
      hist_born[i][ivar]->SetFillColor(kGreen-2);
      hist_box[i][ivar]->SetFillColor(kGreen-1);
      hist_gjet_pf[i][ivar]->SetFillColor(kOrange-2);
      hist_qcd_pf[i][ivar]->SetFillColor(kOrange-3);
      hist_qcd_ff[i][ivar]->SetFillColor(kOrange+2);
      hist_dy[i][ivar]->SetFillColor(38);
      if (!madgraph) {
	hist_gjet_pp[i][ivar]->SetFillColor(kGreen-3);
	hist_qcd_pp[i][ivar]->SetFillColor(kGreen-4);
      }
    }

    hist_bkg_stack[ivar] = new THStack("hist_bkg_stack_"+var[ivar],"Background");
    hist_bkg_stack_sig[ivar] = new THStack("hist_bkg_stack_sig_"+var[ivar],"Background in signal region");
    hist_bkg_stack[ivar]->Add(hist_box[7][ivar]);
    hist_bkg_stack[ivar]->Add(hist_born[7][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_box[3][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_born[3][ivar]);
    if (!madgraph) {
      hist_bkg_stack[ivar]->Add(hist_gjet_pp[7][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_pp[7][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_gjet_pp[3][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_pp[3][ivar]);
    }
    if (fakes) {
      hist_bkg_stack[ivar]->Add(hist_gjet_pf[7][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_pf[7][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_ff[7][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_gjet_pf[3][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_pf[3][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_ff[3][ivar]);
    }    
    hist_bkg_stack[ivar]->Add(hist_dy[7][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_dy[3][ivar]);

    for (int i=0; i<8; i++) {
      hist_bkg[i][ivar] = (TH1*)hist_box[i][ivar]->Clone();
      hist_bkg[i][ivar]->SetFillColor(0);
      hist_bkg[i][ivar]->Add(hist_born[i][ivar]);
      if (!madgraph) {
	hist_bkg[i][ivar]->Add(hist_gjet_pp[i][ivar]);
	hist_bkg[i][ivar]->Add(hist_qcd_pp[i][ivar]);
      }
      if (fakes) {
	hist_bkg[i][ivar]->Add(hist_gjet_pf[i][ivar]);
	hist_bkg[i][ivar]->Add(hist_qcd_pf[i][ivar]);
	hist_bkg[i][ivar]->Add(hist_qcd_ff[i][ivar]);
      }
      hist_bkg[i][ivar]->Add(hist_dy[i][ivar]);
      hist_bkg[i][ivar]->GetXaxis()->SetTitle(title[ivar]);
      hist_bkg[i][ivar]->GetXaxis()->SetTitleSize(0.04);
    }
    float nsig=hist_sig[ivar]->Integral();
    hist_sig_reweight[ivar] = (TH1*)hist_sig[ivar]->Clone();
    hist_sig[ivar]->Scale(hist_bkg[3][ivar]->Integral()/nsig);

    canvas[ivar] = new TCanvas("c_"+var[ivar],var[ivar],1950,1300);
    canvas[ivar]->Divide(3,3);
    canvas[ivar]->SetFillColor(0);


    canvas[ivar]->cd(3);

    TLegend *leg;
    if (ivar==0 || ivar==1 || ivar==4 || ivar==5 || ivar==6 || ivar==7 || (ivar>=10 && ivar!=14)) {
      leg = new TLegend(.55,.55,.87,.87);
    } else if (ivar==14) {
      leg = new TLegend(.32,.1,.68,.35);
    } else {
      leg = new TLegend(.15,.55,.47,.87);
    }
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetTextSize(.035);
    if (data) leg->AddEntry(hist_data[3][ivar],"Data (4.7fb^{-1})");
    leg->AddEntry(hist_sig[ivar],"GluGlu"+mass_str);
    leg->AddEntry(hist_dy[3][ivar],"DYee+Z","F");
    if (fakes) {
      leg->AddEntry(hist_qcd_ff[3][ivar],"QCD fake-fake","F");
      leg->AddEntry(hist_qcd_pf[3][ivar],"QCD prompt-fake","F");
      leg->AddEntry(hist_gjet_pf[3][ivar],"GJet prompt-fake","F");
    }
    if (!madgraph) {
      leg->AddEntry(hist_qcd_pp[3][ivar],"QCD prompt-prompt","F");
      leg->AddEntry(hist_gjet_pp[3][ivar],"GJet prompt-prompt","F");
    }
    if (!madgraph) {
      leg->AddEntry(hist_born[3][ivar],"Pythia Born","F");
    } else {
      leg->AddEntry(hist_born[3][ivar],"Madgraph DiPhotonJets","F");
    }
    leg->AddEntry(hist_box[3][ivar],"Pythia Box","F");


    float max = hist_bkg_stack_sig[ivar]->GetMaximum();
    if (hist_sig[ivar]->GetMaximum()>max) max=hist_sig[ivar]->GetMaximum();
    if (data && hist_data[3][ivar]->GetMaximum()>max) max=hist_data[3][ivar]->GetMaximum();
    hist_bkg_stack_sig[ivar]->SetMaximum(max*1.1);
    if (ivar==15) hist_bkg_stack_sig[ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_bkg_stack_sig[ivar]->SetMaximum(max*1.6);

    hist_bkg_stack_sig[ivar]->Draw("hist");
    if (ivar==1) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    if (ivar==15) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(-1.*sidebandWidth,sidebandWidth);
    //if (ivar==21) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitleSize(0.04);
    if (ivar!=9 && ivar!=3 && ivar!=14) leg->Draw();
    hist_sig[ivar]->Draw("same");
    if (data) hist_data[3][ivar]->Draw("same,e");

    txt = new TLatex();
    txt->SetNDC();
    txt->SetTextSize(0.045);
    if (ivar==2 || ivar==8) {
      txt->DrawLatex(0.15,0.42,"Data in signal region");
      txt->DrawLatex(0.15,0.36,"MC in signal region");
    } else {
      txt->DrawLatex(0.15,0.82,"Data in signal region");
      txt->DrawLatex(0.15,0.76,"MC in signal region");
    }


    canvas[ivar]->cd(1);

    float nbkg_sig = hist_bkg[3][ivar]->Integral();
    hist_bkg_reweight[ivar] = (TH1*)hist_bkg[0][ivar]->Clone();
    hist_bkg_reweight[ivar]->Add(hist_bkg[1][ivar]);
    hist_bkg_reweight[ivar]->Add(hist_bkg[2][ivar]);
    hist_bkg_reweight[ivar]->Add(hist_bkg[4][ivar]);
    hist_bkg_reweight[ivar]->Add(hist_bkg[5][ivar]);
    hist_bkg_reweight[ivar]->Add(hist_bkg[6][ivar]);
    hist_bkg_reweight[ivar]->Scale(nbkg_sig/hist_bkg_reweight[ivar]->Integral());

    int rebinFac = 2;
    if (ivar==15) {
      if (sidebandWidth<0.05) {
	rebinFac = 1;
      } else {
	rebinFac = 5;
      }
    }
    if (ivar==1 || ivar==6 ||ivar==7 || ivar==9) rebinFac=1;
    hist_bkg_reweight[ivar]->Rebin(rebinFac);
    float sf=hist_bkg_reweight[ivar]->Integral();
    hist_bkg_reweight[ivar]->Scale(1./sf);
    for (int i=0; i<8; i++) {
      hist_bkg[i][ivar]->Rebin(rebinFac);
      float sf=hist_bkg[i][ivar]->Integral();
      hist_bkg[i][ivar]->Scale(1./sf);
      if (i!=3) {
	hist_bkg[i][ivar]->SetMarkerStyle(20);
	hist_bkg[i][ivar]->SetMarkerSize(0.8);
      }
    }

    max = hist_bkg[0][ivar]->GetMaximum();
    if (hist_bkg[1][ivar]->GetMaximum()>max) max=hist_bkg[1][ivar]->GetMaximum();
    if (hist_bkg[2][ivar]->GetMaximum()>max) max=hist_bkg[2][ivar]->GetMaximum();
    if (hist_bkg[3][ivar]->GetMaximum()>max) max=hist_bkg[3][ivar]->GetMaximum();
    if (hist_bkg[4][ivar]->GetMaximum()>max) max=hist_bkg[4][ivar]->GetMaximum();
    if (hist_bkg[5][ivar]->GetMaximum()>max) max=hist_bkg[5][ivar]->GetMaximum();
    if (hist_bkg[6][ivar]->GetMaximum()>max) max=hist_bkg[6][ivar]->GetMaximum();
    hist_bkg[3][ivar]->SetMaximum(max*1.1);
    hist_bkg[3][ivar]->SetMinimum(0.);
    if (ivar==15) hist_bkg[3][ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_bkg[3][ivar]->SetMaximum(max*1.6);

    hist_bkg[0][ivar]->SetLineColor(kRed-1);
    hist_bkg[1][ivar]->SetLineColor(kRed-2);
    hist_bkg[2][ivar]->SetLineColor(kRed-3);
    hist_bkg[4][ivar]->SetLineColor(kGreen-3);
    hist_bkg[5][ivar]->SetLineColor(kGreen-2);
    hist_bkg[6][ivar]->SetLineColor(kGreen-1);
    hist_bkg[0][ivar]->SetMarkerColor(kRed-1);
    hist_bkg[1][ivar]->SetMarkerColor(kRed-2);
    hist_bkg[2][ivar]->SetMarkerColor(kRed-3);
    hist_bkg[4][ivar]->SetMarkerColor(kGreen-3);
    hist_bkg[5][ivar]->SetMarkerColor(kGreen-2);
    hist_bkg[6][ivar]->SetMarkerColor(kGreen-1);
    hist_bkg[3][ivar]->SetFillColor(38);
    for (int i=0; i<8; i++) hist_bkg[i][ivar]->SetLineWidth(2);
    hist_bkg[3][ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg[3][ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_bkg[3][ivar]->GetXaxis()->SetRangeUser(0.,7.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg[3][ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    if (ivar==15) hist_bkg[3][ivar]->GetXaxis()->SetRangeUser(-1.*sidebandWidth,sidebandWidth);

    hist_bkg[3][ivar]->Draw("e2");
    if (nSB>2) hist_bkg[0][ivar]->Draw("same");
    if (nSB>1) hist_bkg[1][ivar]->Draw("same");
    hist_bkg[2][ivar]->Draw("same");
    hist_bkg[4][ivar]->Draw("same");
    if (nSB>1) hist_bkg[5][ivar]->Draw("same");
    if (nSB>2) hist_bkg[6][ivar]->Draw("same");

    TLegend *leg2 = (TLegend*)leg->Clone();
    leg2->Clear();
    if (nSB>2) leg2->AddEntry(hist_bkg[0][ivar],"MC Low sideband 3");
    if (nSB>1) leg2->AddEntry(hist_bkg[1][ivar],"MC Low sideband 2");
    leg2->AddEntry(hist_bkg[2][ivar],"MC Low sideband 1");
    leg2->AddEntry(hist_bkg[3][ivar],"MC Signal region");
    leg2->AddEntry(hist_bkg[4][ivar],"MC High sideband 1");
    if (nSB>1) leg2->AddEntry(hist_bkg[5][ivar],"MC High sideband 2");
    if (nSB>2) leg2->AddEntry(hist_bkg[6][ivar],"MC High sideband 3");
    if (ivar!=9) leg2->Draw();


    canvas[ivar]->cd(6);

    hist_data_reweight[ivar] = (TH1*)hist_data[0][ivar]->Clone();
    hist_data_reweight[ivar]->Sumw2();
    hist_data_reweight[ivar]->Add(hist_data[1][ivar]);
    hist_data_reweight[ivar]->Add(hist_data[2][ivar]);
    hist_data_reweight[ivar]->Add(hist_data[4][ivar]);
    hist_data_reweight[ivar]->Add(hist_data[5][ivar]);
    hist_data_reweight[ivar]->Add(hist_data[6][ivar]);
    hist_data_reweight[ivar]->Scale(hist_data[3][ivar]->Integral()/hist_data_reweight[ivar]->Integral());

    //hist_data_reweight[ivar]->Add(hist_data[2][ivar],hist_data[3][ivar]->Integral()/(nSB*2.*hist_data[2][ivar]->Integral())) << endl;
    //hist_data_reweight[ivar]->Add(hist_data[4][ivar],hist_data[3][ivar]->Integral()/(nSB*2.*hist_data[4][ivar]->Integral())) << endl;

    hist_sig_reweight[ivar]->Scale(nbkg_sig/nsig);

    max = hist_bkg_stack_sig[ivar]->GetMaximum();
    if (hist_sig_reweight[ivar]->GetMaximum()>max) max=hist_sig_reweight[ivar]->GetMaximum();
    if (data && hist_data_reweight[ivar]->GetMaximum()>max) max=hist_data_reweight[ivar]->GetMaximum();
    hist_bkg_stack_sig[ivar]->SetMaximum(max*1.1);
    if (ivar==15) hist_bkg_stack_sig[ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_bkg_stack_sig[ivar]->SetMaximum(max*1.6);

    hist_sig_reweight[ivar]->SetLineColor(4);
    hist_sig_reweight[ivar]->SetLineWidth(2.5);
    if (data) hist_data_reweight[ivar]->SetMarkerStyle(20);
    if (data) hist_data_reweight[ivar]->SetMarkerSize(.8);

    hist_bkg_stack_sig[ivar]->Draw("hist");
    if (ivar!=9 && ivar!=3 && ivar!=14) leg->Draw();
    hist_sig_reweight[ivar]->Draw("same");
    if (data) hist_data_reweight[ivar]->Draw("same,e");

    if (ivar==2 || ivar==8) {
      txt->DrawLatex(0.15,0.42,"Background Model (data)");
      txt->DrawLatex(0.15,0.36,"MC in signal region");
    } else {
      txt->DrawLatex(0.15,0.82,"Background Model (data)");
      txt->DrawLatex(0.15,0.76,"MC in signal region");
    }


    canvas[ivar]->cd(4);

    hist_bkg_sig[ivar] = (TH1*)hist_bkg[3][ivar]->Clone();

    max = hist_bkg_reweight[ivar]->GetMaximum();
    if (hist_bkg_sig[ivar]->GetMaximum()>max) max=hist_bkg_sig[ivar]->GetMaximum();
    hist_bkg_reweight[ivar]->SetMaximum(max*1.1);
    hist_bkg_reweight[ivar]->SetMinimum(0.);
    if (ivar==15) hist_bkg_reweight[ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_bkg_reweight[ivar]->SetMaximum(max*1.6);

    hist_bkg_reweight[ivar]->SetMarkerStyle(20);
    hist_bkg_reweight[ivar]->SetMarkerSize(.8);
    hist_bkg_reweight[ivar]->SetMarkerColor(4);
    hist_bkg_reweight[ivar]->SetLineColor(4);
    hist_bkg_reweight[ivar]->SetLineWidth(2);
    hist_bkg_reweight[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg_reweight[ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    if (ivar==15) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(-1.*sidebandWidth,sidebandWidth);
    //if (ivar==21) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    hist_bkg_sig[ivar]->SetFillColor(38);
    hist_bkg_reweight[ivar]->Draw("e");
    hist_bkg_sig[ivar]->Draw("e2,same");
    hist_bkg_reweight[ivar]->Draw("e,same");

    TLegend *leg3 = (TLegend*)leg->Clone();
    leg3->Clear();
    leg3->AddEntry(hist_bkg_sig[ivar],"MC Signal region");
    leg3->AddEntry(hist_bkg_reweight[ivar],"MC Background Model");
    if (ivar!=9) leg3->Draw();


    canvas[ivar]->cd(5);

    hist_data_sig[ivar] = (TH1*)hist_data[3][ivar]->Clone();
    hist_data_sig[ivar]->SetMarkerStyle(0);
    //hist_data_sig[ivar]->Scale(hist_data_reweight[ivar]->Integral()/hist_data_sig[ivar]->Integral());

    hist_data_reweight_rebin[ivar] = (TH1*)hist_data_reweight[ivar]->Clone();
    hist_data_reweight_rebin[ivar]->Rebin(rebinFac);
    hist_data_sig[ivar]->Rebin(rebinFac);

    max = hist_data_reweight_rebin[ivar]->GetMaximum();
    if (hist_data_sig[ivar]->GetMaximum()>max) max=hist_data_sig[ivar]->GetMaximum();
    hist_data_reweight_rebin[ivar]->SetMaximum(max*1.1);
    hist_data_reweight_rebin[ivar]->SetMinimum(0.);
    if (ivar==15) hist_data_reweight_rebin[ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_data_reweight_rebin[ivar]->SetMaximum(max*1.6);

    hist_data_reweight_rebin[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_data_reweight_rebin[ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    if (ivar==15) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(-1.*sidebandWidth,sidebandWidth);
    //if (ivar==21) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    hist_data_sig[ivar]->SetFillColor(38);
    hist_data_reweight_rebin[ivar]->Draw("e");
    hist_data_sig[ivar]->Draw("e2,same");
    hist_data_reweight_rebin[ivar]->Draw("e,same");

    TLegend *leg3 = (TLegend*)leg->Clone();
    leg3->Clear();
    leg3->AddEntry(hist_data_sig[ivar],"Data: Signal region");
    leg3->AddEntry(hist_data_reweight_rebin[ivar],"Data: Background Model");
    if (ivar!=9) leg3->Draw();


    canvas[ivar]->cd(2);

    for (int i=0; i<7; i++) {
      if (i!=3) {
	hist_data[i][ivar]->Rebin(rebinFac);
	float sf=hist_data[i][ivar]->Integral();
	hist_data[i][ivar]->Scale(1./sf);
	hist_data[i][ivar]->SetMarkerStyle(20);
	hist_data[i][ivar]->SetMarkerSize(0.8);
      }
    }

    max = hist_data[0][ivar]->GetMaximum();
    if (hist_data[1][ivar]->GetMaximum()>max) max=hist_data[1][ivar]->GetMaximum();
    if (hist_data[2][ivar]->GetMaximum()>max) max=hist_data[2][ivar]->GetMaximum();
    if (hist_data[4][ivar]->GetMaximum()>max) max=hist_data[4][ivar]->GetMaximum();
    if (hist_data[5][ivar]->GetMaximum()>max) max=hist_data[5][ivar]->GetMaximum();
    if (hist_data[6][ivar]->GetMaximum()>max) max=hist_data[6][ivar]->GetMaximum();
    hist_data[sb_low][ivar]->SetMaximum(max*1.1);
    hist_data[sb_low][ivar]->SetMinimum(0.);
    if (ivar==15) hist_data[sb_low][ivar]->SetMaximum(max*2.);
    if (ivar==6 || ivar==7) hist_data[sb_low][ivar]->SetMaximum(max*1.6);

    hist_data[0][ivar]->SetLineColor(kRed-1);
    hist_data[1][ivar]->SetLineColor(kRed-2);
    hist_data[2][ivar]->SetLineColor(kRed-3);
    hist_data[4][ivar]->SetLineColor(kGreen-3);
    hist_data[5][ivar]->SetLineColor(kGreen-2);
    hist_data[6][ivar]->SetLineColor(kGreen-1);
    hist_data[0][ivar]->SetMarkerColor(kRed-1);
    hist_data[1][ivar]->SetMarkerColor(kRed-2);
    hist_data[2][ivar]->SetMarkerColor(kRed-3);
    hist_data[4][ivar]->SetMarkerColor(kGreen-3);
    hist_data[5][ivar]->SetMarkerColor(kGreen-2);
    hist_data[6][ivar]->SetMarkerColor(kGreen-1);
    for (int i=0; i<8; i++) hist_data[i][ivar]->SetLineWidth(2);
    hist_data[sb_low][ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_data[sb_low][ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_data[sb_low][ivar]->GetXaxis()->SetRangeUser(0.,7.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_data[sb_low][ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    if (ivar==15) hist_data[sb_low][ivar]->GetXaxis()->SetRangeUser(-1.*sidebandWidth,sidebandWidth);

    hist_data[sb_low][ivar]->Draw();
    if (nSB>2) hist_data[1][ivar]->Draw("same");
    if (nSB>1) hist_data[2][ivar]->Draw("same");
    hist_data[4][ivar]->Draw("same");
    if (nSB>1) hist_data[5][ivar]->Draw("same");
    if (nSB>2) hist_data[6][ivar]->Draw("same");

    TLegend *leg2_data = (TLegend*)leg->Clone();
    leg2_data->Clear();
    if (nSB>2) leg2_data->AddEntry(hist_data[0][ivar],"Data Low sideband 3");
    if (nSB>1) leg2_data->AddEntry(hist_data[1][ivar],"Data Low sideband 2");
    leg2_data->AddEntry(hist_data[2][ivar],"Data Low sideband 1");
    leg2_data->AddEntry(hist_data[4][ivar],"Data High sideband 1");
    if (nSB>1) leg2_data->AddEntry(hist_data[5][ivar],"Data High sideband 2");
    if (nSB>2) leg2_data->AddEntry(hist_data[6][ivar],"Data High sideband 3");
    if (ivar!=9) leg2_data->Draw();


    canvas[ivar]->cd(7);

    hist_bkg_ratio[ivar] = (TH1*)hist_bkg_reweight[ivar]->Clone();
    hist_bkg_ratio[ivar]->Reset();
    hist_bkg_ratio[ivar]->Add(hist_bkg_reweight[ivar]);
    hist_bkg_ratio[ivar]->Divide(hist_bkg_sig[ivar]);

    hist_bkg_ratio[ivar]->SetMarkerStyle(20);
    hist_bkg_ratio[ivar]->SetMarkerColor(4);
    hist_bkg_ratio[ivar]->SetMarkerSize(.8);

    hist_bkg_ratio[ivar]->SetMaximum(2.);
    hist_bkg_ratio[ivar]->SetMinimum(0.);
    hist_bkg_ratio[ivar]->Draw();

    float xmin = hist_bkg_ratio[ivar]->GetXaxis()->GetXmin();
    if (ivar==1 || ivar==6 || ivar==7 ||ivar==9) xmin=0.;
    float xmax = hist_bkg_ratio[ivar]->GetXaxis()->GetXmax();
    if (ivar==15) {
      xmin=-1*sidebandWidth;
      xmax=sidebandWidth;
    }
    TLine *line2 = new TLine(xmin,1.,xmax,1.);
    line2->SetLineColor(4);
    line2->SetLineWidth(2);
    line2->Draw();

    txt->DrawLatex(0.15,0.82, "MC Background Model / MC in signal region");


    canvas[ivar]->cd(8);
    
    hist_data_ratio[ivar] = (TH1*)hist_data_reweight_rebin[ivar]->Clone();
    hist_data_ratio[ivar]->Reset();
    hist_data_ratio[ivar]->Add(hist_data_reweight_rebin[ivar]);
    hist_data_ratio[ivar]->Divide(hist_data_sig[ivar]);

    hist_data_ratio[ivar]->SetMarkerStyle(20);
    hist_data_ratio[ivar]->SetMarkerSize(.8);

    hist_data_ratio[ivar]->SetMaximum(2.);
    hist_data_ratio[ivar]->SetMinimum(0.);
    hist_data_ratio[ivar]->Draw("e");

    line2->Draw();

    txt->DrawLatex(0.15,0.82,"Background Model / Data in signal region");

    canvas[ivar]->SaveAs(outdir+var[ivar]+".gif");

  }

  //if (mass_in != 110) {
    for (int j=22; j<24; j++) {

      hist_sig[j]->SetLineColor(2);
      hist_sig[j]->SetLineWidth(2.5);

      for (int i=0; i<7; i++) {
	hist_data[i][j]->SetMarkerStyle(20);
	hist_data[i][j]->SetMarkerSize(.8);
	hist_data[i][j]->SetLineWidth(2);
	hist_data[i][j]->GetXaxis()->SetTitle(title[j]);
	hist_data[i][j]->GetXaxis()->SetTitleSize(0.04);
	hist_bkg[i][j]->SetLineWidth(2);
	hist_bkg[i][j]->GetXaxis()->SetTitle(title[j]);
	hist_bkg[i][j]->GetXaxis()->SetTitleSize(0.04);
	if (i!=3) {
	  hist_bkg[i][j]->SetMarkerStyle(20);
	  hist_bkg[i][j]->SetMarkerSize(0.8);
	}
      }

      hist_bkg[0][j]->SetLineColor(kRed-1);
      hist_bkg[1][j]->SetLineColor(kRed-2);
      hist_bkg[2][j]->SetLineColor(kRed-3);
      hist_bkg[3][j]->SetLineColor(4);
      hist_bkg[4][j]->SetLineColor(kGreen-3);
      hist_bkg[5][j]->SetLineColor(kGreen-2);
      hist_bkg[6][j]->SetLineColor(kGreen-1);
      hist_bkg[0][j]->SetMarkerColor(kRed-1);
      hist_bkg[1][j]->SetMarkerColor(kRed-2);
      hist_bkg[2][j]->SetMarkerColor(kRed-3);
      hist_bkg[3][j]->SetMarkerColor(4);
      hist_bkg[4][j]->SetMarkerColor(kGreen-3);
      hist_bkg[5][j]->SetMarkerColor(kGreen-2);
      hist_bkg[6][j]->SetMarkerColor(kGreen-1);
      hist_bkg[3][j]->SetFillColor(38);
      hist_bkgModel[j]->SetMarkerStyle(20);
      hist_bkgModel[j]->SetMarkerSize(.8);
      hist_bkgModel[j]->SetLineColor(4);
      hist_bkgModel[j]->SetFillColor(38);
      hist_bkgModel[j]->SetLineWidth(2);
      hist_bkgModel[j]->GetXaxis()->SetTitle(title[j]);
      hist_bkgModel[j]->GetXaxis()->SetTitleSize(0.04);
      hist_sig_x5 = (TH1*)hist_sig[j]->Clone();
      hist_sig_x10 = (TH1*)hist_sig[j]->Clone();
      hist_sig_x5->SetLineStyle(2);
      hist_sig[j]->SetLineStyle(3);
      hist_sig[j]->Scale(1.);
      hist_sig_x5->Scale(5.);
      hist_sig_x10->Scale(10.);

      nbkg_sig = hist_bkg[3][j]->Integral();
      for (int i=0; i<7; i++) {
	float sf=nbkg_sig/hist_bkg[i][j]->Integral();
	hist_bkg_scaled[i][j] = (TH1*)hist_bkg[i][j]->Clone();
	hist_bkg_scaled[i][j]->Scale(sf);

	sf=Ndata_sig/hist_data[i][j]->Integral();
	hist_data_scaled[i][j] = (TH1*)hist_data[i][j]->Clone();
	hist_data_scaled[i][j]->Scale(sf);

      }

      float nbkg_sig = hist_bkg[3][j]->Integral();
      hist_bkg_reweight[j] = (TH1*)hist_bkg[0][j]->Clone();
      hist_bkg_reweight[j]->Add(hist_bkg[1][j]);
      hist_bkg_reweight[j]->Add(hist_bkg[2][j]);
      hist_bkg_reweight[j]->Add(hist_bkg[4][j]);
      hist_bkg_reweight[j]->Add(hist_bkg[5][j]);
      hist_bkg_reweight[j]->Add(hist_bkg[6][j]);
      hist_bkg_reweight[j]->Scale(nbkg_sig/hist_bkg_reweight[j]->Integral());

      /*
      hist_data_reweight[j] = (TH1*)hist_data[0][j]->Clone();
      hist_data_reweight[j]->Add(hist_data[1][j]);
      hist_data_reweight[j]->Add(hist_data[2][j]);
      hist_data_reweight[j]->Add(hist_data[4][j]);
      hist_data_reweight[j]->Add(hist_data[5][j]);
      hist_data_reweight[j]->Add(hist_data[6][j]);
      hist_data_reweight[j]->Scale(Ndata_sig/hist_data_reweight[j]->Integral());
      hist_data_reweight[j]->SetLineColor(3);
      hist_data_reweight[j]->SetMarkerColor(3);
      */

      canvas[j] = new TCanvas("c_"+var[j],var[j],1950,1300);
      canvas[j]->Divide(3,3);
      canvas[j]->SetFillColor(0);

      canvas[j]->cd(1);
      if (logy) gPad->SetLogy();

      float max = hist_bkgModel[j]->GetMaximum();
      if (hist_data[3][j]->GetMaximum()>max) max=hist_data[3][j]->GetMaximum();
      if (!sob) max=nperbin*2.4/1.1;
      hist_bkgModel[j]->SetMaximum(max*1.1);
      hist_bkgModel[j]->SetMinimum(0.);

      hist_bkgModel[j]->Draw("hist");
      if (logy) hist_bkgModel[j]->GetYaxis()->UnZoom();
      //hist_data_reweight[j]->Draw("same");
      hist_sig[j]->Draw("hist,same");
      hist_sig_x5->Draw("hist,same");
      hist_sig_x10->Draw("hist,same");
      hist_data[3][j]->Draw("same,e");

      TLegend *leg1 = new TLegend(.55,.55,.87,.87);
      leg1->SetBorderSize(0);
      leg1->SetFillColor(10);
      leg1->SetTextSize(.035);
      leg1->AddEntry(hist_data[3][j],"Data (4.7fb^{-1})");
      leg1->AddEntry(hist_sig_x10,"Signal ("+mass_str+" GeV) x1, x5, x10");
      leg1->AddEntry(hist_bkgModel[j],"Background Model","F");
      leg1->Draw();
      leg1->Draw();

      canvas[j]->cd(5);

      hist_data_highMinusLow = (TH1*)hist_data_scaled[sb_high][j]->Clone();
      //hist_data_highMinusLow->Add(hist_data_scaled[sb_low][j],-1.);
      hist_data_highMinusLow->Divide(hist_data_scaled[sb_low][j]);
      hist_data_highMinusLow->SetMarkerColor(4);
      hist_data_highMinusLow->SetLineColor(4);
      hist_data_highMinusLow->SetMaximum(2.);
      hist_data_highMinusLow->SetMinimum(.0);
      hist_data_highMinusLow->Draw("e");

      float xmin = hist_data[1][j]->GetXaxis()->GetXmin();
      float xmax = hist_data[1][j]->GetXaxis()->GetXmax();
      TLine *line1 = new TLine(xmin,1.,xmax,1.);
      line1->SetLineColor(4);
      line1->SetLineWidth(2);
      line1->Draw();

      txt->DrawLatex(0.15,0.82, "highest sideband / lowest sideband");

      TF1 *f1 = new TF1("f1","[0]",0.,xmax);
      hist_data_highMinusLow_fit = (TH1*)hist_data_highMinusLow->Clone();
      hist_data_highMinusLow_fit->Fit(f1,"Q0");
      Float_t chi2 = hist_data_highMinusLow_fit->GetFunction("f1")->GetChisquare()/hist_data_highMinusLow_fit->GetFunction("f1")->GetNDF();
      Float_t chi2_rounded = floor(chi2*100.+0.5)/100.;
      TString chi2_str = "#chi^{2} = ";
      chi2_str+=chi2;
      gStyle->SetPaintTextFormat("3.2g");
      cout << chi2_str << endl;
      //txt->DrawLatex(0.7,0.2, chi2_str);

      /*
	canvas[j]->cd(3);

	hist_bkg[3][j]->SetMaximum(max*1.1);
	hist_bkg[3][j]->SetMinimum(0.);

	hist_bkg[3][j]->Draw("hist");
	hist_sig[j]->Draw("hist,same");
	hist_data[3][j]->Draw("same,e");

	leg3 = (TLegend*)leg1->Clone();
	leg3->Clear();
	leg3->AddEntry(hist_data[3][j],"Data (4.7fb^{-1})");
	leg3->AddEntry(hist_sig[j],"Signal ("+mass_str+" GeV) x 10");
	leg3->AddEntry(hist_bkg[3][j],"Background MC","F");
	leg3->Draw();
      */

      canvas[j]->cd(6);

      hist_bkg_highMinusLow = (TH1*)hist_bkg_scaled[sb_high][j]->Clone();
      //hist_bkg_highMinusLow->Add(hist_bkg_scaled[sb_low][j],-1.);
      hist_bkg_highMinusLow->Divide(hist_bkg_scaled[sb_low][j]);
      hist_bkg_highMinusLow->SetMarkerColor(4);
      hist_bkg_highMinusLow->SetLineColor(4);
      hist_bkg_highMinusLow->SetMaximum(2.);
      hist_bkg_highMinusLow->SetMinimum(0.);
      hist_bkg_highMinusLow->Draw("e");

      txt->DrawLatex(0.15,0.82, "MC highest sideband / MC lowest sideband");

      line1->Draw();

      canvas[j]->cd(4);

      hist_dataMinusModel = (TH1*)hist_data[3][j]->Clone();
      hist_dataMinusModel->Add(hist_bkgModel[j],-1.);

      float max1=100.;
      float min1=-100.;
      if (hist_dataMinusModel->GetMaximum()>max1) max1=hist_dataMinusModel->GetMaximum()*1.2;
      if (hist_dataMinusModel->GetMinimum()<min1) min1=hist_dataMinusModel->GetMinimum()*1.2;
      if (max1>fabs(min1)) min1=-1.*max1;
      if (fabs(min1)>max1) max1=-1.*min1;
      //if (hist_sig_x10[j]->GetMaximum()>max1) max1=hist_sig_x10[j]->GetMaximum()*1.05;
      hist_dataMinusModel->SetMaximum(max1);
      hist_dataMinusModel->SetMinimum(min1);
      hist_dataMinusModel->Draw("e");
      hist_sig[j]->Draw("hist,same");
      hist_sig_x5->Draw("hist,same");
      hist_sig_x10->Draw("hist,same");

      TLegend *leg2 = new TLegend(.45,.15,.87,.32);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(10);
      leg2->SetTextSize(.035);
      leg2->AddEntry(hist_dataMinusModel,"Data (4.7fb^{-1}) - background model","LP");
      leg2->AddEntry(hist_sig_x10,"Signal ("+mass_str+" GeV) x1, x5, x10");
      leg2->Draw();

      TLine *line1a = new TLine(xmin,0.,xmax,0.);
      line1a->SetLineColor(4);
      line1a->SetLineWidth(2);
      line1a->Draw();


      canvas[j]->cd(7);

      hist_SoverB = (TH1*)hist_sig[j]->Clone();
      hist_SoverB->Divide(hist_bkgModel[j]);
      //hist_SoverB->Divide(hist_balg[j]);
      hist_SoverB->GetYaxis()->SetTitle("S/B");
      hist_SoverB->GetXaxis()->SetTitle(title[j]);
      hist_SoverB->SetMarkerStyle(20);
      hist_SoverB->SetMarkerSize(1.);
      hist_SoverB->SetMarkerColor(2);
      hist_SoverB->SetLineStyle(1);
      hist_SoverB->SetLineWidth(2.);
      hist_SoverB->Draw();


      canvas[j]->cd(2);
      if (logy) gPad->SetLogy();

      hist_data_scaled[0][j]->SetLineColor(kRed-1);
      hist_data_scaled[1][j]->SetLineColor(kRed-2);
      hist_data_scaled[2][j]->SetLineColor(kRed-3);
      hist_data_scaled[4][j]->SetLineColor(kGreen-3);
      hist_data_scaled[5][j]->SetLineColor(kGreen-2);
      hist_data_scaled[6][j]->SetLineColor(kGreen-1);
      hist_data_scaled[0][j]->SetMarkerColor(kRed-1);
      hist_data_scaled[1][j]->SetMarkerColor(kRed-2);
      hist_data_scaled[2][j]->SetMarkerColor(kRed-3);
      hist_data_scaled[4][j]->SetMarkerColor(kGreen-3);
      hist_data_scaled[5][j]->SetMarkerColor(kGreen-2);
      hist_data_scaled[6][j]->SetMarkerColor(kGreen-1);

      hist_data_scaled[sb_low][j]->SetMaximum(max*1.1);
      hist_data_scaled[sb_low][j]->SetMinimum(0.);
      hist_data_scaled[sb_low][j]->Draw("e");
      if (logy) hist_data_scaled[sb_low][j]->GetYaxis()->UnZoom();
      if (nSB>2) hist_data_scaled[1][j]->Draw("e,same");
      if (nSB>1) hist_data_scaled[2][j]->Draw("e,same");
      hist_data_scaled[4][j]->Draw("e,same");
      if (nSB>1) hist_data_scaled[5][j]->Draw("e,same");
      if (nSB>2) hist_data_scaled[6][j]->Draw("e,same");

      leg6 = (TLegend*)leg1->Clone();
      leg6->Clear();
      if (nSB>2) leg6->AddEntry(hist_data_scaled[6][j],"High sideband 3","LP");
      if (nSB>1) leg6->AddEntry(hist_data_scaled[5][j],"High sideband 2","LP");
      leg6->AddEntry(hist_data_scaled[4][j],"High sideband 1","LP");
      leg6->AddEntry(hist_data_scaled[2][j],"Low sideband 1","LP");
      if (nSB>1) leg6->AddEntry(hist_data_scaled[1][j],"Low sideband 2","LP");
      if (nSB>2) leg6->AddEntry(hist_data_scaled[0][j],"Low sideband 3","LP");
      leg6->Draw();

      canvas[j]->cd(9);
      if (logy) gPad->SetLogy();

      hist_bkg_reweight[j]->SetMaximum(max*1.1);
      hist_bkg_reweight[j]->SetMinimum(0.);
      hist_bkg_reweight[j]->SetLineColor(4);
      hist_bkg_reweight[j]->SetMarkerColor(4);
      hist_bkg[3][j]->SetFillColor(38);

      hist_bkg_reweight[j]->Draw("e");
      if (logy) hist_bkg_reweight[j]->GetYaxis()->UnZoom();
      hist_bkg[3][j]->Draw("e2,same");
      hist_bkg_reweight[j]->Draw("e,same");

      TLegend *leg7 = (TLegend*)leg1->Clone();
      leg7->Clear();
      leg7->AddEntry(hist_bkg[3][j],"MC Signal region","F");
      leg7->AddEntry(hist_bkg_reweight[j],"MC Background Model");
      leg7->Draw();

      canvas[j]->cd(3);
      if (logy) gPad->SetLogy();

      hist_bkg_scaled[3][j]->SetMaximum(max*1.1);
      hist_bkg_scaled[3][j]->SetMinimum(0.);

      hist_bkg_scaled[3][j]->Draw("e2");
      if (logy) hist_bkg_scaled[3][j]->GetYaxis()->UnZoom();
      if (nSB>2) hist_bkg_scaled[0][j]->Draw("e,same");
      if (nSB>1) hist_bkg_scaled[1][j]->Draw("e,same");
      hist_bkg_scaled[2][j]->Draw("e,same");
      hist_bkg_scaled[4][j]->Draw("e,same");
      if (nSB>1) hist_bkg_scaled[5][j]->Draw("e,same");
      if (nSB>2) hist_bkg_scaled[6][j]->Draw("e,same");

      leg8 = (TLegend*)leg1->Clone();
      leg8->Clear();
      if (nSB>2) leg8->AddEntry(hist_bkg_scaled[6][j],"MC High sideband 3","LP");
      if (nSB>1) leg8->AddEntry(hist_bkg_scaled[5][j],"MC High sideband 2","LP");
      leg8->AddEntry(hist_bkg_scaled[4][j],"MC High sideband 1","LP");
      leg8->AddEntry(hist_bkg_scaled[3][j],"MC Signal region","F");
      leg8->AddEntry(hist_bkg_scaled[2][j],"MC Low sideband 1","LP");
      if (nSB>1) leg8->AddEntry(hist_bkg_scaled[1][j],"MC Low sideband 2","LP");
      if (nSB>2) leg8->AddEntry(hist_bkg_scaled[0][j],"MC Low sideband 3","LP");
      leg8->Draw();

      if (!rebin) {
	canvas[j]->SaveAs(outdir+var[j]+sob_str+".gif");
      } else {
	canvas[j]->SaveAs(outdir+var[j]+sob_str+"_nominalbins.gif");
      }

    }
    //}

/*
  hist_mass_born->Scale(bornSF);
  hist_mass_box->Scale(1.3);
  hist_mass_gjet_pf->Scale(1.3);
  hist_mass_qcd_pf->Scale(1.3);
  hist_mass_qcd_ff->Scale(1.);
  hist_mass_dy->Scale(1.15*2321./992.);
  if (!madgraph) {
    hist_mass_gjet_pp->Scale(1.3);
    hist_mass_qcd_pp->Scale(1.3);
  }
*/

  hist_mass_sig->SetLineColor(4);
  hist_mass_sig->SetLineWidth(2.5);
  if (data) hist_mass_data->SetMarkerStyle(20);
  if (data) hist_mass_data->SetMarkerSize(.8);
  hist_mass_born->SetFillColor(kGreen-2);
  hist_mass_box->SetFillColor(kGreen-1);
  hist_mass_gjet_pf->SetFillColor(kOrange-2);
  hist_mass_qcd_pf->SetFillColor(kOrange-3);
  hist_mass_qcd_ff->SetFillColor(kOrange+2);
  hist_mass_dy->SetFillColor(38);
  if (!madgraph) {
    hist_mass_gjet_pp->SetFillColor(kGreen-3);
    hist_mass_qcd_pp->SetFillColor(kGreen-4);
  }

  hist_mass_bkg_stack = new THStack("hist_mass_bkg_stack","Background");
  hist_mass_bkg_stack->Add(hist_mass_box);
  hist_mass_bkg_stack->Add(hist_mass_born);
  if (!madgraph) {
    hist_mass_bkg_stack->Add(hist_mass_gjet_pp);
    hist_mass_bkg_stack->Add(hist_mass_qcd_pp);
  }
  if (fakes) {
    hist_mass_bkg_stack->Add(hist_mass_gjet_pf);
    hist_mass_bkg_stack->Add(hist_mass_qcd_pf);
    hist_mass_bkg_stack->Add(hist_mass_qcd_ff);
  }    
  hist_mass_bkg_stack->Add(hist_mass_dy);
  
  hist_mass_bkg = (TH1*)hist_mass_box->Clone();
  hist_mass_bkg->Add(hist_mass_born);
  if (!madgraph) {
    hist_mass_bkg->Add(hist_mass_gjet_pp);
    hist_mass_bkg->Add(hist_mass_qcd_pp);
  }
  if (fakes) {
    hist_mass_bkg->Add(hist_mass_gjet_pf);
    hist_mass_bkg->Add(hist_mass_qcd_pf);
    hist_mass_bkg->Add(hist_mass_qcd_ff);
  }
  hist_mass_bkg->Add(hist_mass_dy);
  hist_mass_sig->Scale(10.);

  cout << "Data/MC scale factor = " << hist_mass_data->Integral(hist_mass_data->FindBin(sideband_boundaries[sb_low]),hist_mass_data->FindBin(sideband_boundaries[sb_high]))/hist_mass_bkg->Integral(hist_mass_data->FindBin(sideband_boundaries[sb_low]),hist_mass_data->FindBin(sideband_boundaries[sb_high])) << endl;

  TCanvas *c_mgg = new TCanvas("c_mgg","Mgg",1000,700);
  c_mgg->SetFillColor(0);

  float max = hist_mass_bkg_stack->GetMaximum();
  if (hist_mass_sig->GetMaximum()>max) max=hist_mass_sig->GetMaximum();
  if (data && hist_mass_data->GetMaximum()>max) max=hist_mass_data->GetMaximum();
  //max=380.;
  hist_mass_bkg_stack->SetMaximum(max*1.05);

  hist_mass_bkg_stack->Draw();
  hist_mass_bkg_stack->GetXaxis()->SetRangeUser(90.,190.);

  TBox* sideband_box[2];
  sideband_box[0] = new TBox(sideband_boundaries[sb_low],0.,sideband_boundaries[3],max*1.1);
  sideband_box[1] = new TBox(sideband_boundaries[4],0.,sideband_boundaries[sb_high+1],max*1.1);
  for (int i=0; i<2; i++) {
    sideband_box[i]->SetFillColor(38);
    sideband_box[i]->SetFillStyle(3002);
    sideband_box[i]->Draw("same");
  }

  TBox* signalRegion_box = new TBox(signalRegion_boundaries[0],0.,signalRegion_boundaries[1],max*1.1);
  signalRegion_box->SetFillColor(kRed);
  signalRegion_box->SetFillStyle(3002);
  signalRegion_box->Draw("same");

  hist_mass_bkg_stack->Draw("same");
  //hist_mass_bkg_stack->GetXaxis()->SetRangeUser(95.,185.);
  hist_mass_bkg_stack->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  hist_mass_bkg_stack->GetXaxis()->SetTitleSize(0.04);
  hist_mass_sig->Draw("same");
  if (data) hist_mass_data->Draw("same,e");
  gPad->RedrawAxis();

  TLine* line_sig[2];
  for (int i=0; i<2; i++) {
    line_sig[i] = new TLine(signalRegion_boundaries[i],0.,signalRegion_boundaries[i],max*1.1);
    line_sig[i]->SetLineColor(kRed-7);
    line_sig[i]->SetLineWidth(3);
    line_sig[i]->SetLineStyle(9);
    line_sig[i]->Draw();
  }

  TLine* line[8];
  for (int i=0; i<8; i++) {
    line[i] = new TLine(sideband_boundaries[i],0.,sideband_boundaries[i],max*1.1);
    line[i]->SetLineColor(38);
    line[i]->SetLineWidth(3);
    line[i]->SetLineStyle(9);
  }
  if (nSB>2) line[0]->Draw();
  if (nSB>1) line[1]->Draw("same");
  line[2]->Draw();
  line[3]->Draw();
  line[4]->Draw();
  line[5]->Draw();
  if (nSB>1) line[6]->Draw("same");
  if (nSB>2) line[7]->Draw("same");

  TLegend *leg_mass = new TLegend(.65,.6,.87,.87);
  leg_mass->SetBorderSize(0);
  leg_mass->SetFillColor(10);
  leg_mass->SetTextSize(.025);
  if (data) leg_mass->AddEntry(hist_mass_data,"Data (4.7fb^{-1})");
  leg_mass->AddEntry(hist_mass_sig,"GluGlu"+mass_str+" #times10");
  leg_mass->AddEntry(hist_mass_dy,"DYee+Z","F");
  if (fakes) {
    leg_mass->AddEntry(hist_mass_qcd_ff,"QCD fake-fake","F");
    leg_mass->AddEntry(hist_mass_qcd_pf,"QCD prompt-fake","F");
    leg_mass->AddEntry(hist_mass_gjet_pf,"GJet prompt-fake","F");
  }
  if (!madgraph) {
    leg_mass->AddEntry(hist_mass_qcd_pp,"QCD prompt-prompt","F");
    leg_mass->AddEntry(hist_mass_gjet_pp,"GJet prompt-prompt","F");
  }
  if (!madgraph) {
    leg_mass->AddEntry(hist_mass_born,"Pythia Born","F");
  } else {
    leg_mass->AddEntry(hist_mass_born,"Madgraph DiPhotonJets","F");
  }
  leg_mass->AddEntry(hist_mass_box,"Pythia Box","F");
  leg_mass->Draw();

  c_mgg->SaveAs(outdir+"mass.gif");

  /*
  TCanvas *c_res = new TCanvas("c_res","Mass Resolution",1800,1800);
  c_res->SetFillColor(0);
  c_res->Divide(3,4);

  for (int icat=0; icat<9; icat++) {

    if (icat==0) c_res->cd(1);
    if (icat==1) c_res->cd(2);
    if (icat==2) c_res->cd(3);
    if (icat==3) c_res->cd(5);
    if (icat==4) c_res->cd(6);
    if (icat==5) c_res->cd(8);
    if (icat==6) c_res->cd(9);
    if (icat==7) c_res->cd(11);
    if (icat==8) c_res->cd(12);

//     hist_res_born[icat]->Scale(bornSF*dataMC_sf);
//     hist_res_box[icat]->Scale(1.3*dataMC_sf);
//     hist_res_gjet_pf[icat]->Scale(1.3*dataMC_sf);
//     hist_res_qcd_pf[icat]->Scale(1.3*dataMC_sf);
//     hist_res_qcd_ff[icat]->Scale(1.*dataMC_sf);
//     hist_res_dy[icat]->Scale(1.15*2321./992.*dataMC_sf);
//     if (!madgraph) {
//       hist_res_gjet_pp[icat]->Scale(1.3*dataMC_sf);
//       hist_res_qcd_pp[icat]->Scale(1.3*dataMC_sf);
//     }

    hist_res_born[icat]->Scale(dataMC_sf);
    hist_res_box[icat]->Scale(dataMC_sf);
    hist_res_gjet_pf[icat]->Scale(dataMC_sf);
    hist_res_qcd_pf[icat]->Scale(dataMC_sf);
    hist_res_qcd_ff[icat]->Scale(dataMC_sf);
    hist_res_dy[icat]->Scale(dataMC_sf);
    if (!madgraph) {
      hist_res_gjet_pp[icat]->Scale(dataMC_sf);
      hist_res_qcd_pp[icat]->Scale(dataMC_sf);
    }

    hist_res_sig[icat]->SetLineColor(4);
    hist_res_sig[icat]->SetLineWidth(2.5);
    if (data) hist_res_data[icat]->SetMarkerStyle(20);
    if (data) hist_res_data[icat]->SetMarkerSize(.5);
    hist_res_born[icat]->SetFillColor(kGreen-2);
    hist_res_box[icat]->SetFillColor(kGreen-1);
    hist_res_gjet_pf[icat]->SetFillColor(kOrange-2);
    hist_res_qcd_pf[icat]->SetFillColor(kOrange-3);
    hist_res_qcd_ff[icat]->SetFillColor(kOrange+2);
    hist_res_dy[icat]->SetFillColor(38);
    if (!madgraph) {
      hist_res_gjet_pp[icat]->SetFillColor(kGreen-3);
      hist_res_qcd_pp[icat]->SetFillColor(kGreen-4);
    }

    hist_res_bkg_stack[icat] = new THStack("hist_res_bkg_stack","Background");
    hist_res_bkg_stack[icat]->Add(hist_res_box[icat]);
    hist_res_bkg_stack[icat]->Add(hist_res_born[icat]);
    if (!madgraph) {
      hist_res_bkg_stack[icat]->Add(hist_res_gjet_pp[icat]);
      hist_res_bkg_stack[icat]->Add(hist_res_qcd_pp[icat]);
    }
    if (fakes) {
      hist_res_bkg_stack[icat]->Add(hist_res_gjet_pf[icat]);
      hist_res_bkg_stack[icat]->Add(hist_res_qcd_pf[icat]);
      hist_res_bkg_stack[icat]->Add(hist_res_qcd_ff[icat]);
    }    
    hist_res_bkg_stack[icat]->Add(hist_res_dy[icat]);
  
    hist_res_bkg[icat] = (TH1*)hist_res_box[icat]->Clone();
    hist_res_bkg[icat]->Add(hist_res_born[icat]);
    if (!madgraph) {
      hist_res_bkg[icat]->Add(hist_res_gjet_pp[icat]);
      hist_res_bkg[icat]->Add(hist_res_qcd_pp[icat]);
    }
    if (fakes) {
      hist_res_bkg[icat]->Add(hist_res_gjet_pf[icat]);
      hist_res_bkg[icat]->Add(hist_res_qcd_pf[icat]);
      hist_res_bkg[icat]->Add(hist_res_qcd_ff[icat]);
    }
    hist_res_bkg[icat]->Add(hist_res_dy[icat]);
    hist_res_sig[icat]->Scale(hist_res_bkg[icat]->Integral()/hist_res_sig[icat]->Integral());

    float max = hist_res_bkg_stack[icat]->GetMaximum();
    if (hist_res_sig[icat]->GetMaximum()>max) max=hist_res_sig[icat]->GetMaximum();
    if (data && hist_res_data[icat]->GetMaximum()>max) max=hist_res_data[icat]->GetMaximum();
    hist_res_bkg_stack[icat]->SetMaximum(max*1.05);

    hist_res_bkg_stack[icat]->Draw();
    hist_res_bkg_stack[icat]->GetXaxis()->SetRangeUser(0.,10.);
    hist_res_bkg_stack[icat]->GetXaxis()->SetTitle("#sigma_{M}");
    hist_res_bkg_stack[icat]->GetXaxis()->SetTitleSize(0.04);
    hist_res_sig[icat]->Draw("same");
    if (data) hist_res_data[icat]->Draw("same,e");

    txt_cat = new TLatex();
    txt_cat->SetNDC();
    txt_cat->SetTextSize(0.06);
    if (icat==0) {
      txt->DrawLatex(0.4,0.8,"All");
    } else {
      if (icat<5) {
	txt->DrawLatex(0.4,0.8,"High p_{T}^{#gamma#gamma}");
      } else {
	txt->DrawLatex(0.4,0.8,"Low p_{T}^{#gamma#gamma}");
      }
      if (icat==1 || icat==2 || icat==5 || icat==6) {
	txt->DrawLatex(0.4,0.75,"Both EB");
      } else {
	txt->DrawLatex(0.4,0.75,"One EE");
      }
      if (icat%2==1) {
	txt->DrawLatex(0.4,0.7,"High R9");
      } else {
	txt->DrawLatex(0.4,0.7,"Low R9");
      }
    }

    leg_mass->Draw();
  }

  c_res->SaveAs(outdir+"sigmaM_cat.gif");
  */
}

void SetHistogramErrors(TH1* h) {

  if (h->GetEntries()==0 || h->Integral()==0) return;

  float sf = h->Integral()/h->GetEntries();
  for (int ibin=0; ibin<h->GetNbinsX()+1; ibin++) {
    float err = sqrt(h->GetBinContent(ibin)/sf);
    h->SetBinError(ibin,err*sf);
  }

}

std::vector<double> optimizedReverseBinning(TH1F *hb,int nTargetBins,bool revise_target, bool use_n_target){
    // Return a set of bins which are "smoother" 

    if (revise_target) {
        if (use_n_target){
           std::cerr << "WARNING -- RooContainer::OptimizedBinning -- Can't use number of Entries as target in revised binning algo " << std::endl; 
           use_n_target = false;  // geometric algo always use revised number of bins, not number of entries
        
        }
    }

    int nBins = hb->GetNbinsX();
    std::vector<double> binEdges;

    double targetNumbers;
    if (use_n_target) targetNumbers = nTargetBins; 
    else targetNumbers = hb->Integral()/nTargetBins;

    if (hb->Integral() < 2*targetNumbers){
        std::cout << "RooContainer::OptimizedBinning -- Not enough entries in histogram for target numbers calculated - " 
              << targetNumbers 
              << ", Returning current bin boundaries "  << std::endl;
        //for (int j=2;j<=nBins+1;j++) binEdges.push_back(hb->GetBinLowEdge(j));
        binEdges.push_back(hb->GetBinLowEdge(1));
        binEdges.push_back(hb->GetBinLowEdge(nBins+1));
        return binEdges;
    }
    binEdges.push_back(hb->GetBinLowEdge(nBins+1));

    std::cout << "RooContainer::optimizedBinning -- Performing Reverse Optimize Binning" <<std::endl;
    double sumBin = 0;
    int i=nBins;
    while (i>=1){
        if (revise_target) targetNumbers = hb->Integral(1,i)/nTargetBins;
        sumBin=hb->GetBinContent(i);
        double highEdge=hb->GetBinLowEdge(i);

        bool carryOn = sumBin <= targetNumbers;
        while ( carryOn){
            if (i>1){
              sumBin+=hb->GetBinContent(i-1);
              highEdge = hb->GetBinLowEdge(i-1);
              carryOn =(sumBin <targetNumbers && i>=1);
              i--;
            } else {
              highEdge = hb->GetBinLowEdge(i);
              carryOn=false;
            }
        }
            binEdges.push_back(highEdge);
        i--;
    }
        if (sumBin < 10) binEdges.erase(binEdges.end()-2);
    reverse(binEdges.begin(),binEdges.end());
    return binEdges;

}
