#include "TMath.h"
#include <iomanip>

/*
This macro can be executed directly in root to produce a set of plots for 
mH=120.  A script to generate and execute macros for all mass hypotheses is 
runBackgroundModelPlots.sh.  This requires that the script 
make_bkgplot_html.sh is run first in order to create a directory structure 
and html pages to organize the plots that are produced.
*/

void backgroundModelPlots(bool www=false, TString outdirname="BDTplots_all", bool madgraph=true, bool fakes=true, bool data=true) {

  int mass_in=120;
  bool equalBinWidths=1;
  bool rebinBdtOut=0;

  const float signalRegionWidth=0.07;
  const float dataMC_sf = 1.04;

  TString mass_str;
  mass_str+=mass_in;
  TString outdir;
  if (www) {
    outdir = "/afs/cern.ch/user/f/futyand/www/hgg/"+outdirname+"/"+mass_str+"/gifs/";
  } else {
    outdir = outdirname+"/";
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  TFile *f_sig = TFile::Open("CMS-HGG_signal_all.root");
  if (data) TFile *f_data = TFile::Open("CMS-HGG_data.root");
  TFile *f_born;
  TFile *f_gjet_pp;
  TFile *f_qcd_pp;
  if (!madgraph) {
    f_born = TFile::Open("CMS-HGG_Born25.root");
    f_gjet_pp = TFile::Open("CMS-HGG_GJetPP.root");
    f_qcd_pp = TFile::Open("CMS-HGG_QCDPP.root");
  }else {
    f_born = TFile::Open("CMS-HGG_DiPhotonJets.root");
  }
  TFile *f_box = TFile::Open("CMS-HGG_Box25.root");
  TFile *f_gjet_pf = TFile::Open("CMS-HGG_GJetPF.root");
  TFile *f_qcd_pf = TFile::Open("CMS-HGG_QCDPF.root");
  TFile *f_qcd_ff = TFile::Open("CMS-HGG_QCDFF.root");
  TFile *f_dy = TFile::Open("CMS-HGG_DyeeMG.root");

  TFile *f_bdtout = TFile::Open("CMS-HGG_1658pb_mva_20-09-11.root");

  TH1* hist_sig[23];
  TH1* hist_data[4][23];
  TH1* hist_born[4][23];
  TH1* hist_box[4][23];
  TH1* hist_gjet_pp[4][23];
  TH1* hist_gjet_pf[4][23];
  TH1* hist_qcd_pp[4][23];
  TH1* hist_qcd_pf[4][23];
  TH1* hist_qcd_ff[4][23];
  TH1* hist_dy[4][23];
  TH1* hist_bkg[4][23];
  TH1* hist_bkgModel[23];
  TH1* hist_bkg_scaled[4][23];
  TH1* hist_data_scaled[4][23];
  TH1* hist_bkg_sig[23];
  TH1* hist_data_sig[23];
  TH1* hist_sig_reweight[23];
  TH1* hist_bkg_reweight[23];
  TH1* hist_data_reweight[23];
  TH1* hist_data_reweight_rebin[23];
  TH1* hist_bkg_ratio[23];
  TH1* hist_bkg_ratio_low[23];
  TH1* hist_bkg_ratio_high[23];
  TH1* hist_data_ratio[23];
  THStack* hist_bkg_stack[23];
  THStack* hist_bkg_stack_sig[23];

  TH2 *hist2D_sig;
  TH2 *hist2D_data;
  TH2 *hist2D_born;
  TH2 *hist2D_box;
  TH2 *hist2D_gjet_pp;
  TH2 *hist2D_gjet_pf;
  TH2 *hist2D_qcd_pp;
  TH2 *hist2D_qcd_pf;
  TH2 *hist2D_qcd_ff;
  TH2 *hist2D_bkg;

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

  TH1 *hist_res_sig[9];
  TH1 *hist_res_data[9];
  TH1 *hist_res_born[9];
  TH1 *hist_res_box[9];
  TH1 *hist_res_gjet_pp[9];
  TH1 *hist_res_gjet_pf[9];
  TH1 *hist_res_qcd_pp[9];
  TH1 *hist_res_qcd_pf[9];
  TH1 *hist_res_qcd_ff[9];
  TH1 *hist_res_dy[9];
  TH1 *hist_res_bkg[9];
  THStack* hist_res_bkg_stack[9];

  f_bdtout->cd();

  if (!equalBinWidths) {

    hist_sig[22] = (TH1*)th1f_sig_grad_120_cat0->Clone();
    hist_data[1][22] = (TH1*)th1f_bkg_low_grad_120_cat0->Clone();
    hist_data[2][22] = (TH1*)th1f_data_grad_120_cat0->Clone();
    hist_data[3][22] = (TH1*)th1f_bkg_high_grad_120_cat0->Clone();
    hist_bkgModel[22] = (TH1*)th1f_bkg_grad_120_cat0->Clone();
    hist_bkg[1][22] = (TH1*)th1f_bkg_mc_low_grad_120_cat0->Clone();
    hist_bkg[2][22] = (TH1*)th1f_bkg_mc_grad_120_cat0->Clone();
    hist_bkg[3][22] = (TH1*)th1f_bkg_mc_high_grad_120_cat0->Clone();
    hist_sig[22]->GetYaxis()->SetTitle("");
    hist_bkgModel[22]->GetYaxis()->SetTitle("");
    for (int i=1; i<4; i++) {
      hist_data[i][22]->GetYaxis()->SetTitle("");
      hist_bkg[i][22]->GetYaxis()->SetTitle("");
    }

  } else {

    int nbins=th1f_sig_grad_120_cat0->GetNbinsX();
    if (rebinBdtOut) nbins-=1;
    hist_sig[22] = new TH1F("hist_sig","hist_sig",nbins,0,float(nbins));
    hist_data[1][22] = new TH1F("hist_data_low","hist_data_low",nbins,0,float(nbins));
    hist_data[2][22] = new TH1F("hist_data_sig","hist_data_sig",nbins,0,float(nbins));
    hist_data[3][22] = new TH1F("hist_data_high","hist_data_high",nbins,0,float(nbins));
    hist_bkgModel[22] = new TH1F("hist_bkgModel","hist_bkgModel",nbins,0,float(nbins));
    hist_bkg[1][22] = new TH1F("hist_bkg_low","hist_bkg_low",nbins,0,float(nbins));
    hist_bkg[2][22] = new TH1F("hist_bkg_sig","hist_bkg_sig",nbins,0,float(nbins));
    hist_bkg[3][22] = new TH1F("hist_bkg_high","hist_bkg_high",nbins,0,float(nbins));
    for (int ibin=0; ibin<nbins+1; ibin++) {
      hist_sig[22]->SetBinContent(ibin,th1f_sig_grad_120_cat0->GetBinContent(ibin));
      hist_sig[22]->SetBinError(ibin,th1f_sig_grad_120_cat0->GetBinError(ibin));
      hist_data[1][22]->SetBinContent(ibin,th1f_bkg_low_grad_120_cat0->GetBinContent(ibin));
      hist_data[1][22]->SetBinError(ibin,th1f_bkg_low_grad_120_cat0->GetBinError(ibin));
      hist_data[2][22]->SetBinContent(ibin,th1f_data_grad_120_cat0->GetBinContent(ibin));
      hist_data[2][22]->SetBinError(ibin,th1f_data_grad_120_cat0->GetBinError(ibin));
      hist_data[3][22]->SetBinContent(ibin,th1f_bkg_high_grad_120_cat0->GetBinContent(ibin));
      hist_data[3][22]->SetBinError(ibin,th1f_bkg_high_grad_120_cat0->GetBinError(ibin));
      hist_bkgModel[22]->SetBinError(ibin,th1f_bkg_grad_120_cat0->GetBinError(ibin));
      hist_bkgModel[22]->SetBinContent(ibin,th1f_bkg_grad_120_cat0->GetBinContent(ibin));
      hist_bkg[1][22]->SetBinContent(ibin,th1f_bkg_mc_low_grad_120_cat0->GetBinContent(ibin));
      hist_bkg[1][22]->SetBinError(ibin,th1f_bkg_mc_low_grad_120_cat0->GetBinError(ibin));
      hist_bkg[2][22]->SetBinContent(ibin,th1f_bkg_mc_grad_120_cat0->GetBinContent(ibin));
      hist_bkg[2][22]->SetBinError(ibin,th1f_bkg_mc_grad_120_cat0->GetBinError(ibin));
      hist_bkg[3][22]->SetBinContent(ibin,th1f_bkg_mc_high_grad_120_cat0->GetBinContent(ibin));
      hist_bkg[3][22]->SetBinError(ibin,th1f_bkg_mc_high_grad_120_cat0->GetBinError(ibin));
    }

    if (rebinBdtOut) {
      hist_sig[22]->Rebin(2);
      for (int i=1; i<4; i++) {
	hist_data[i][22]->Rebin(2);
	hist_bkg[i][22]->Rebin(2);
	hist_bkgModel[22]->Rebin(2);
      }
    }

  }

  f_sig->cd();

  hist_sig[0] = (TH1*)ptOverM_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[1] = (TH1*)eta_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[2] = (TH1*)deltaPhi_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[3] = (TH1*)helicityAngle_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[4] = (TH1*)pho1_pt_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[5] = (TH1*)pho2_pt_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[6] = (TH1*)pho1_eta_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[7] = (TH1*)pho2_eta_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[8] = (TH1*)pho_minr9_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[9] = (TH1*)maxeta_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[10] = (TH1*)ptOverM_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[11] = (TH1*)pho1_ptOverM_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[12] = (TH1*)pho2_ptOverM_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[13] = (TH1*)ht_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[14] = (TH1*)deltaEta_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[15] = (TH1*)deltaMOverMH_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[16] = (TH1*)pho1_ptOverMH_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[17] = (TH1*)pho2_ptOverMH_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[18] = (TH1*)sigmaMOverM_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[20] = (TH1*)deltaMOverSigmaM_msig_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_sig[21] = (TH1*)sigmaM_msig_cat2_gluglu_H_gg_120_pu2011->Clone();

  hist2D_sig = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_mass_sig = (TH1*)all_mass_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[0] = (TH1*)sigmaM_cat0_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[1] = (TH1*)sigmaM_cat1_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[2] = (TH1*)sigmaM_cat2_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[3] = (TH1*)sigmaM_cat3_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[4] = (TH1*)sigmaM_cat4_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[5] = (TH1*)sigmaM_cat5_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[6] = (TH1*)sigmaM_cat6_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[7] = (TH1*)sigmaM_cat7_gluglu_H_gg_120_pu2011->Clone();
  hist_res_sig[8] = (TH1*)sigmaM_cat8_gluglu_H_gg_120_pu2011->Clone();

  if (data) {

    f_data->cd();
    hist_data[1][0] = (TH1*)ptOverM_mlow_cat2_Data->Clone();
    hist_data[1][1] = (TH1*)eta_mlow_cat2_Data->Clone();
    hist_data[1][2] = (TH1*)deltaPhi_mlow_cat2_Data->Clone();
    hist_data[1][3] = (TH1*)helicityAngle_mlow_cat2_Data->Clone();
    hist_data[1][4] = (TH1*)pho1_pt_mlow_cat2_Data->Clone();
    hist_data[1][5] = (TH1*)pho2_pt_mlow_cat2_Data->Clone();
    hist_data[1][6] = (TH1*)pho1_eta_mlow_cat2_Data->Clone();
    hist_data[1][7] = (TH1*)pho2_eta_mlow_cat2_Data->Clone();
    hist_data[1][8] = (TH1*)pho_minr9_mlow_cat2_Data->Clone();
    hist_data[1][9] = (TH1*)maxeta_mlow_cat2_Data->Clone();
    hist_data[1][10] = (TH1*)ptOverM_mlow_cat2_Data->Clone();
    hist_data[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_Data->Clone();
    hist_data[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_Data->Clone();
    hist_data[1][13] = (TH1*)ht_mlow_cat2_Data->Clone();
    hist_data[1][14] = (TH1*)deltaEta_mlow_cat2_Data->Clone();
    hist_data[1][15] = (TH1*)deltaMOverMH_mlow_cat2_Data->Clone();
    hist_data[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_Data->Clone();
    hist_data[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_Data->Clone();
    hist_data[1][18] = (TH1*)sigmaMOverM_mlow_cat2_Data->Clone();
    hist_data[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_Data->Clone();
    hist_data[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_Data->Clone();
    hist_data[1][21] = (TH1*)sigmaM_mlow_cat2_Data->Clone();

    hist_data[2][0] = (TH1*)ptOverM_msig_cat2_Data->Clone();
    hist_data[2][1] = (TH1*)eta_msig_cat2_Data->Clone();
    hist_data[2][2] = (TH1*)deltaPhi_msig_cat2_Data->Clone();
    hist_data[2][3] = (TH1*)helicityAngle_msig_cat2_Data->Clone();
    hist_data[2][4] = (TH1*)pho1_pt_msig_cat2_Data->Clone();
    hist_data[2][5] = (TH1*)pho2_pt_msig_cat2_Data->Clone();
    hist_data[2][6] = (TH1*)pho1_eta_msig_cat2_Data->Clone();
    hist_data[2][7] = (TH1*)pho2_eta_msig_cat2_Data->Clone();
    hist_data[2][8] = (TH1*)pho_minr9_msig_cat2_Data->Clone();
    hist_data[2][9] = (TH1*)maxeta_msig_cat2_Data->Clone();
    hist_data[2][10] = (TH1*)ptOverM_msig_cat2_Data->Clone();
    hist_data[2][11] = (TH1*)pho1_ptOverM_msig_cat2_Data->Clone();
    hist_data[2][12] = (TH1*)pho2_ptOverM_msig_cat2_Data->Clone();
    hist_data[2][13] = (TH1*)ht_msig_cat2_Data->Clone();
    hist_data[2][14] = (TH1*)deltaEta_msig_cat2_Data->Clone();
    hist_data[2][15] = (TH1*)deltaMOverMH_msig_cat2_Data->Clone();
    hist_data[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_Data->Clone();
    hist_data[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_Data->Clone();
    hist_data[2][18] = (TH1*)sigmaMOverM_msig_cat2_Data->Clone();
    hist_data[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_Data->Clone();
    hist_data[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_Data->Clone();
    hist_data[2][21] = (TH1*)sigmaM_msig_cat2_Data->Clone();

    hist_data[3][0] = (TH1*)ptOverM_mhigh_cat2_Data->Clone();
    hist_data[3][1] = (TH1*)eta_mhigh_cat2_Data->Clone();
    hist_data[3][2] = (TH1*)deltaPhi_mhigh_cat2_Data->Clone();
    hist_data[3][3] = (TH1*)helicityAngle_mhigh_cat2_Data->Clone();
    hist_data[3][4] = (TH1*)pho1_pt_mhigh_cat2_Data->Clone();
    hist_data[3][5] = (TH1*)pho2_pt_mhigh_cat2_Data->Clone();
    hist_data[3][6] = (TH1*)pho1_eta_mhigh_cat2_Data->Clone();
    hist_data[3][7] = (TH1*)pho2_eta_mhigh_cat2_Data->Clone();
    hist_data[3][8] = (TH1*)pho_minr9_mhigh_cat2_Data->Clone();
    hist_data[3][9] = (TH1*)maxeta_mhigh_cat2_Data->Clone();
    hist_data[3][10] = (TH1*)ptOverM_mhigh_cat2_Data->Clone();
    hist_data[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_Data->Clone();
    hist_data[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_Data->Clone();
    hist_data[3][13] = (TH1*)ht_mhigh_cat2_Data->Clone();
    hist_data[3][14] = (TH1*)deltaEta_mhigh_cat2_Data->Clone();
    hist_data[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_Data->Clone();
    hist_data[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_Data->Clone();
    hist_data[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_Data->Clone();
    hist_data[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_Data->Clone();
    hist_data[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_Data->Clone();
    hist_data[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_Data->Clone();
    hist_data[3][21] = (TH1*)sigmaM_mhigh_cat2_Data->Clone();

    hist2D_data = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_Data->Clone();
    hist_mass_data = (TH1*)all_mass_cat0_Data->Clone();
    hist_res_data[0] = (TH1*)sigmaM_cat0_Data->Clone();
    hist_res_data[1] = (TH1*)sigmaM_cat1_Data->Clone();
    hist_res_data[2] = (TH1*)sigmaM_cat2_Data->Clone();
    hist_res_data[3] = (TH1*)sigmaM_cat3_Data->Clone();
    hist_res_data[4] = (TH1*)sigmaM_cat4_Data->Clone();
    hist_res_data[5] = (TH1*)sigmaM_cat5_Data->Clone();
    hist_res_data[6] = (TH1*)sigmaM_cat6_Data->Clone();
    hist_res_data[7] = (TH1*)sigmaM_cat7_Data->Clone();
    hist_res_data[8] = (TH1*)sigmaM_cat8_Data->Clone();

  }

  if (!madgraph) {

    f_born->cd();
    hist_born[1][0] = (TH1*)ptOverM_mlow_cat2_Born25->Clone();
    hist_born[1][1] = (TH1*)eta_mlow_cat2_Born25->Clone();
    hist_born[1][2] = (TH1*)deltaPhi_mlow_cat2_Born25->Clone();
    hist_born[1][3] = (TH1*)helicityAngle_mlow_cat2_Born25->Clone();
    hist_born[1][4] = (TH1*)pho1_pt_mlow_cat2_Born25->Clone();
    hist_born[1][5] = (TH1*)pho2_pt_mlow_cat2_Born25->Clone();
    hist_born[1][6] = (TH1*)pho1_eta_mlow_cat2_Born25->Clone();
    hist_born[1][7] = (TH1*)pho2_eta_mlow_cat2_Born25->Clone();
    hist_born[1][8] = (TH1*)pho_minr9_mlow_cat2_Born25->Clone();
    hist_born[1][9] = (TH1*)maxeta_mlow_cat2_Born25->Clone();
    hist_born[1][10] = (TH1*)ptOverM_mlow_cat2_Born25->Clone();
    hist_born[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_Born25->Clone();
    hist_born[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_Born25->Clone();
    hist_born[1][13] = (TH1*)ht_mlow_cat2_Born25->Clone();
    hist_born[1][14] = (TH1*)deltaEta_mlow_cat2_Born25->Clone();
    hist_born[1][15] = (TH1*)deltaMOverMH_mlow_cat2_Born25->Clone();
    hist_born[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_Born25->Clone();
    hist_born[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_Born25->Clone();
    hist_born[1][18] = (TH1*)sigmaMOverM_mlow_cat2_Born25->Clone();
    hist_born[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_Born25->Clone();
    hist_born[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_Born25->Clone();
    hist_born[1][21] = (TH1*)sigmaM_mlow_cat2_Born25->Clone();

    hist_born[2][0] = (TH1*)ptOverM_msig_cat2_Born25->Clone();
    hist_born[2][1] = (TH1*)eta_msig_cat2_Born25->Clone();
    hist_born[2][2] = (TH1*)deltaPhi_msig_cat2_Born25->Clone();
    hist_born[2][3] = (TH1*)helicityAngle_msig_cat2_Born25->Clone();
    hist_born[2][4] = (TH1*)pho1_pt_msig_cat2_Born25->Clone();
    hist_born[2][5] = (TH1*)pho2_pt_msig_cat2_Born25->Clone();
    hist_born[2][6] = (TH1*)pho1_eta_msig_cat2_Born25->Clone();
    hist_born[2][7] = (TH1*)pho2_eta_msig_cat2_Born25->Clone();
    hist_born[2][8] = (TH1*)pho_minr9_msig_cat2_Born25->Clone();
    hist_born[2][9] = (TH1*)maxeta_msig_cat2_Born25->Clone();
    hist_born[2][10] = (TH1*)ptOverM_msig_cat2_Born25->Clone();
    hist_born[2][11] = (TH1*)pho1_ptOverM_msig_cat2_Born25->Clone();
    hist_born[2][12] = (TH1*)pho2_ptOverM_msig_cat2_Born25->Clone();
    hist_born[2][13] = (TH1*)ht_msig_cat2_Born25->Clone();
    hist_born[2][14] = (TH1*)deltaEta_msig_cat2_Born25->Clone();
    hist_born[2][15] = (TH1*)deltaMOverMH_msig_cat2_Born25->Clone();
    hist_born[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_Born25->Clone();
    hist_born[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_Born25->Clone();
    hist_born[2][18] = (TH1*)sigmaMOverM_msig_cat2_Born25->Clone();
    hist_born[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_Born25->Clone();
    hist_born[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_Born25->Clone();
    hist_born[2][21] = (TH1*)sigmaM_msig_cat2_Born25->Clone();

    hist_born[3][0] = (TH1*)ptOverM_mhigh_cat2_Born25->Clone();
    hist_born[3][1] = (TH1*)eta_mhigh_cat2_Born25->Clone();
    hist_born[3][2] = (TH1*)deltaPhi_mhigh_cat2_Born25->Clone();
    hist_born[3][3] = (TH1*)helicityAngle_mhigh_cat2_Born25->Clone();
    hist_born[3][4] = (TH1*)pho1_pt_mhigh_cat2_Born25->Clone();
    hist_born[3][5] = (TH1*)pho2_pt_mhigh_cat2_Born25->Clone();
    hist_born[3][6] = (TH1*)pho1_eta_mhigh_cat2_Born25->Clone();
    hist_born[3][7] = (TH1*)pho2_eta_mhigh_cat2_Born25->Clone();
    hist_born[3][8] = (TH1*)pho_minr9_mhigh_cat2_Born25->Clone();
    hist_born[3][9] = (TH1*)maxeta_mhigh_cat2_Born25->Clone();
    hist_born[3][10] = (TH1*)ptOverM_mhigh_cat2_Born25->Clone();
    hist_born[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_Born25->Clone();
    hist_born[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_Born25->Clone();
    hist_born[3][13] = (TH1*)ht_mhigh_cat2_Born25->Clone();
    hist_born[3][14] = (TH1*)deltaEta_mhigh_cat2_Born25->Clone();
    hist_born[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_Born25->Clone();
    hist_born[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_Born25->Clone();
    hist_born[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_Born25->Clone();
    hist_born[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_Born25->Clone();
    hist_born[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_Born25->Clone();
    hist_born[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_Born25->Clone();
    hist_born[3][21] = (TH1*)sigmaM_mhigh_cat2_Born25->Clone();

    hist2D_born = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_Born25->Clone();
    hist_mass_born = (TH1*)all_mass_cat0_Born25->Clone();
    hist_res_born[0] = (TH1*)sigmaM_cat0_Born25->Clone();
    hist_res_born[1] = (TH1*)sigmaM_cat1_Born25->Clone();
    hist_res_born[2] = (TH1*)sigmaM_cat2_Born25->Clone();
    hist_res_born[3] = (TH1*)sigmaM_cat3_Born25->Clone();
    hist_res_born[4] = (TH1*)sigmaM_cat4_Born25->Clone();
    hist_res_born[5] = (TH1*)sigmaM_cat5_Born25->Clone();
    hist_res_born[6] = (TH1*)sigmaM_cat6_Born25->Clone();
    hist_res_born[7] = (TH1*)sigmaM_cat7_Born25->Clone();
    hist_res_born[8] = (TH1*)sigmaM_cat8_Born25->Clone();

  } else {

    f_born->cd();
    hist_born[1][0] = (TH1*)ptOverM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][1] = (TH1*)eta_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][2] = (TH1*)deltaPhi_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][3] = (TH1*)helicityAngle_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][4] = (TH1*)pho1_pt_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][5] = (TH1*)pho2_pt_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][6] = (TH1*)pho1_eta_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][7] = (TH1*)pho2_eta_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][8] = (TH1*)pho_minr9_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][9] = (TH1*)maxeta_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][10] = (TH1*)ptOverM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][13] = (TH1*)ht_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][14] = (TH1*)deltaEta_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][15] = (TH1*)deltaMOverMH_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][18] = (TH1*)sigmaMOverM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_DiPhotonJets->Clone();
    hist_born[1][21] = (TH1*)sigmaM_mlow_cat2_DiPhotonJets->Clone();

    hist_born[2][0] = (TH1*)ptOverM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][1] = (TH1*)eta_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][2] = (TH1*)deltaPhi_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][3] = (TH1*)helicityAngle_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][4] = (TH1*)pho1_pt_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][5] = (TH1*)pho2_pt_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][6] = (TH1*)pho1_eta_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][7] = (TH1*)pho2_eta_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][8] = (TH1*)pho_minr9_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][9] = (TH1*)maxeta_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][10] = (TH1*)ptOverM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][11] = (TH1*)pho1_ptOverM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][12] = (TH1*)pho2_ptOverM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][13] = (TH1*)ht_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][14] = (TH1*)deltaEta_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][15] = (TH1*)deltaMOverMH_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][18] = (TH1*)sigmaMOverM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_DiPhotonJets->Clone();
    hist_born[2][21] = (TH1*)sigmaM_msig_cat2_DiPhotonJets->Clone();

    hist_born[3][0] = (TH1*)ptOverM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][1] = (TH1*)eta_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][2] = (TH1*)deltaPhi_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][3] = (TH1*)helicityAngle_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][4] = (TH1*)pho1_pt_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][5] = (TH1*)pho2_pt_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][6] = (TH1*)pho1_eta_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][7] = (TH1*)pho2_eta_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][8] = (TH1*)pho_minr9_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][9] = (TH1*)maxeta_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][10] = (TH1*)ptOverM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][13] = (TH1*)ht_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][14] = (TH1*)deltaEta_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_DiPhotonJets->Clone();
    hist_born[3][21] = (TH1*)sigmaM_mhigh_cat2_DiPhotonJets->Clone();

    hist2D_born = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_DiPhotonJets->Clone();
    hist_mass_born = (TH1*)all_mass_cat0_DiPhotonJets->Clone();
    hist_res_born[0] = (TH1*)sigmaM_cat0_DiPhotonJets->Clone();
    hist_res_born[1] = (TH1*)sigmaM_cat1_DiPhotonJets->Clone();
    hist_res_born[2] = (TH1*)sigmaM_cat2_DiPhotonJets->Clone();
    hist_res_born[3] = (TH1*)sigmaM_cat3_DiPhotonJets->Clone();
    hist_res_born[4] = (TH1*)sigmaM_cat4_DiPhotonJets->Clone();
    hist_res_born[5] = (TH1*)sigmaM_cat5_DiPhotonJets->Clone();
    hist_res_born[6] = (TH1*)sigmaM_cat6_DiPhotonJets->Clone();
    hist_res_born[7] = (TH1*)sigmaM_cat7_DiPhotonJets->Clone();
    hist_res_born[8] = (TH1*)sigmaM_cat8_DiPhotonJets->Clone();
  }

  f_box->cd();
  hist_box[1][0] = (TH1*)ptOverM_mlow_cat2_Box25->Clone();
  hist_box[1][1] = (TH1*)eta_mlow_cat2_Box25->Clone();
  hist_box[1][2] = (TH1*)deltaPhi_mlow_cat2_Box25->Clone();
  hist_box[1][3] = (TH1*)helicityAngle_mlow_cat2_Box25->Clone();
  hist_box[1][4] = (TH1*)pho1_pt_mlow_cat2_Box25->Clone();
  hist_box[1][5] = (TH1*)pho2_pt_mlow_cat2_Box25->Clone();
  hist_box[1][6] = (TH1*)pho1_eta_mlow_cat2_Box25->Clone();
  hist_box[1][7] = (TH1*)pho2_eta_mlow_cat2_Box25->Clone();
  hist_box[1][8] = (TH1*)pho_minr9_mlow_cat2_Box25->Clone();
  hist_box[1][9] = (TH1*)maxeta_mlow_cat2_Box25->Clone();
  hist_box[1][10] = (TH1*)ptOverM_mlow_cat2_Box25->Clone();
  hist_box[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_Box25->Clone();
  hist_box[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_Box25->Clone();
  hist_box[1][13] = (TH1*)ht_mlow_cat2_Box25->Clone();
  hist_box[1][14] = (TH1*)deltaEta_mlow_cat2_Box25->Clone();
  hist_box[1][15] = (TH1*)deltaMOverMH_mlow_cat2_Box25->Clone();
  hist_box[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_Box25->Clone();
  hist_box[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_Box25->Clone();
  hist_box[1][18] = (TH1*)sigmaMOverM_mlow_cat2_Box25->Clone();
  hist_box[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_Box25->Clone();
  hist_box[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_Box25->Clone();
  hist_box[1][21] = (TH1*)sigmaM_mlow_cat2_Box25->Clone();

  hist_box[2][0] = (TH1*)ptOverM_msig_cat2_Box25->Clone();
  hist_box[2][1] = (TH1*)eta_msig_cat2_Box25->Clone();
  hist_box[2][2] = (TH1*)deltaPhi_msig_cat2_Box25->Clone();
  hist_box[2][3] = (TH1*)helicityAngle_msig_cat2_Box25->Clone();
  hist_box[2][4] = (TH1*)pho1_pt_msig_cat2_Box25->Clone();
  hist_box[2][5] = (TH1*)pho2_pt_msig_cat2_Box25->Clone();
  hist_box[2][6] = (TH1*)pho1_eta_msig_cat2_Box25->Clone();
  hist_box[2][7] = (TH1*)pho2_eta_msig_cat2_Box25->Clone();
  hist_box[2][8] = (TH1*)pho_minr9_msig_cat2_Box25->Clone();
  hist_box[2][9] = (TH1*)maxeta_msig_cat2_Box25->Clone();
  hist_box[2][10] = (TH1*)ptOverM_msig_cat2_Box25->Clone();
  hist_box[2][11] = (TH1*)pho1_ptOverM_msig_cat2_Box25->Clone();
  hist_box[2][12] = (TH1*)pho2_ptOverM_msig_cat2_Box25->Clone();
  hist_box[2][13] = (TH1*)ht_msig_cat2_Box25->Clone();
  hist_box[2][14] = (TH1*)deltaEta_msig_cat2_Box25->Clone();
  hist_box[2][15] = (TH1*)deltaMOverMH_msig_cat2_Box25->Clone();
  hist_box[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_Box25->Clone();
  hist_box[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_Box25->Clone();
  hist_box[2][18] = (TH1*)sigmaMOverM_msig_cat2_Box25->Clone();
  hist_box[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_Box25->Clone();
  hist_box[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_Box25->Clone();
  hist_box[2][21] = (TH1*)sigmaM_msig_cat2_Box25->Clone();

  hist_box[3][0] = (TH1*)ptOverM_mhigh_cat2_Box25->Clone();
  hist_box[3][1] = (TH1*)eta_mhigh_cat2_Box25->Clone();
  hist_box[3][2] = (TH1*)deltaPhi_mhigh_cat2_Box25->Clone();
  hist_box[3][3] = (TH1*)helicityAngle_mhigh_cat2_Box25->Clone();
  hist_box[3][4] = (TH1*)pho1_pt_mhigh_cat2_Box25->Clone();
  hist_box[3][5] = (TH1*)pho2_pt_mhigh_cat2_Box25->Clone();
  hist_box[3][6] = (TH1*)pho1_eta_mhigh_cat2_Box25->Clone();
  hist_box[3][7] = (TH1*)pho2_eta_mhigh_cat2_Box25->Clone();
  hist_box[3][8] = (TH1*)pho_minr9_mhigh_cat2_Box25->Clone();
  hist_box[3][9] = (TH1*)maxeta_mhigh_cat2_Box25->Clone();
  hist_box[3][10] = (TH1*)ptOverM_mhigh_cat2_Box25->Clone();
  hist_box[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_Box25->Clone();
  hist_box[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_Box25->Clone();
  hist_box[3][13] = (TH1*)ht_mhigh_cat2_Box25->Clone();
  hist_box[3][14] = (TH1*)deltaEta_mhigh_cat2_Box25->Clone();
  hist_box[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_Box25->Clone();
  hist_box[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_Box25->Clone();
  hist_box[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_Box25->Clone();
  hist_box[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_Box25->Clone();
  hist_box[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_Box25->Clone();
  hist_box[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_Box25->Clone();
  hist_box[3][21] = (TH1*)sigmaM_mhigh_cat2_Box25->Clone();

  hist2D_box = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_Box25->Clone();
  hist_mass_box = (TH1*)all_mass_cat0_Box25->Clone();
  hist_res_box[0] = (TH1*)sigmaM_cat0_Box25->Clone();
  hist_res_box[1] = (TH1*)sigmaM_cat1_Box25->Clone();
  hist_res_box[2] = (TH1*)sigmaM_cat2_Box25->Clone();
  hist_res_box[3] = (TH1*)sigmaM_cat3_Box25->Clone();
  hist_res_box[4] = (TH1*)sigmaM_cat4_Box25->Clone();
  hist_res_box[5] = (TH1*)sigmaM_cat5_Box25->Clone();
  hist_res_box[6] = (TH1*)sigmaM_cat6_Box25->Clone();
  hist_res_box[7] = (TH1*)sigmaM_cat7_Box25->Clone();
  hist_res_box[8] = (TH1*)sigmaM_cat8_Box25->Clone();

  if (!madgraph) {

    f_gjet_pp->cd();
    hist_gjet_pp[1][0] = (TH1*)ptOverM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][1] = (TH1*)eta_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][2] = (TH1*)deltaPhi_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][3] = (TH1*)helicityAngle_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][4] = (TH1*)pho1_pt_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][5] = (TH1*)pho2_pt_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][6] = (TH1*)pho1_eta_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][7] = (TH1*)pho2_eta_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][8] = (TH1*)pho_minr9_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][9] = (TH1*)maxeta_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][10] = (TH1*)ptOverM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][13] = (TH1*)ht_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][14] = (TH1*)deltaEta_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][15] = (TH1*)deltaMOverMH_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][18] = (TH1*)sigmaMOverM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_GJetPP->Clone();
    hist_gjet_pp[1][21] = (TH1*)sigmaM_mlow_cat2_GJetPP->Clone();

    hist_gjet_pp[2][0] = (TH1*)ptOverM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][1] = (TH1*)eta_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][2] = (TH1*)deltaPhi_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][3] = (TH1*)helicityAngle_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][4] = (TH1*)pho1_pt_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][5] = (TH1*)pho2_pt_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][6] = (TH1*)pho1_eta_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][7] = (TH1*)pho2_eta_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][8] = (TH1*)pho_minr9_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][9] = (TH1*)maxeta_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][10] = (TH1*)ptOverM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][11] = (TH1*)pho1_ptOverM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][12] = (TH1*)pho2_ptOverM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][13] = (TH1*)ht_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][14] = (TH1*)deltaEta_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][15] = (TH1*)deltaMOverMH_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][18] = (TH1*)sigmaMOverM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_GJetPP->Clone();
    hist_gjet_pp[2][21] = (TH1*)sigmaM_msig_cat2_GJetPP->Clone();

    hist_gjet_pp[3][0] = (TH1*)ptOverM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][1] = (TH1*)eta_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][2] = (TH1*)deltaPhi_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][3] = (TH1*)helicityAngle_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][4] = (TH1*)pho1_pt_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][5] = (TH1*)pho2_pt_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][6] = (TH1*)pho1_eta_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][7] = (TH1*)pho2_eta_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][8] = (TH1*)pho_minr9_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][9] = (TH1*)maxeta_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][10] = (TH1*)ptOverM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][13] = (TH1*)ht_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][14] = (TH1*)deltaEta_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_GJetPP->Clone();
    hist_gjet_pp[3][21] = (TH1*)sigmaM_mhigh_cat2_GJetPP->Clone();

    hist2D_gjet_pp = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_GJetPP->Clone();
    hist_mass_gjet_pp = (TH1*)all_mass_cat0_GJetPP->Clone();
    hist_res_gjet_pp[0] = (TH1*)sigmaM_cat0_GJetPP->Clone();
    hist_res_gjet_pp[1] = (TH1*)sigmaM_cat1_GJetPP->Clone();
    hist_res_gjet_pp[2] = (TH1*)sigmaM_cat2_GJetPP->Clone();
    hist_res_gjet_pp[3] = (TH1*)sigmaM_cat3_GJetPP->Clone();
    hist_res_gjet_pp[4] = (TH1*)sigmaM_cat4_GJetPP->Clone();
    hist_res_gjet_pp[5] = (TH1*)sigmaM_cat5_GJetPP->Clone();
    hist_res_gjet_pp[6] = (TH1*)sigmaM_cat6_GJetPP->Clone();
    hist_res_gjet_pp[7] = (TH1*)sigmaM_cat7_GJetPP->Clone();
    hist_res_gjet_pp[8] = (TH1*)sigmaM_cat8_GJetPP->Clone();

  }

  f_gjet_pf->cd();
  hist_gjet_pf[1][0] = (TH1*)ptOverM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][1] = (TH1*)eta_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][2] = (TH1*)deltaPhi_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][3] = (TH1*)helicityAngle_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][4] = (TH1*)pho1_pt_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][5] = (TH1*)pho2_pt_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][6] = (TH1*)pho1_eta_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][7] = (TH1*)pho2_eta_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][8] = (TH1*)pho_minr9_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][9] = (TH1*)maxeta_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][10] = (TH1*)ptOverM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][13] = (TH1*)ht_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][14] = (TH1*)deltaEta_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][15] = (TH1*)deltaMOverMH_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][18] = (TH1*)sigmaMOverM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_GJetPF->Clone();
  hist_gjet_pf[1][21] = (TH1*)sigmaM_mlow_cat2_GJetPF->Clone();

  hist_gjet_pf[2][0] = (TH1*)ptOverM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][1] = (TH1*)eta_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][2] = (TH1*)deltaPhi_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][3] = (TH1*)helicityAngle_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][4] = (TH1*)pho1_pt_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][5] = (TH1*)pho2_pt_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][6] = (TH1*)pho1_eta_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][7] = (TH1*)pho2_eta_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][8] = (TH1*)pho_minr9_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][9] = (TH1*)maxeta_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][10] = (TH1*)ptOverM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][11] = (TH1*)pho1_ptOverM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][12] = (TH1*)pho2_ptOverM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][13] = (TH1*)ht_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][14] = (TH1*)deltaEta_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][15] = (TH1*)deltaMOverMH_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][18] = (TH1*)sigmaMOverM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_GJetPF->Clone();
  hist_gjet_pf[2][21] = (TH1*)sigmaM_msig_cat2_GJetPF->Clone();

  hist_gjet_pf[3][0] = (TH1*)ptOverM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][1] = (TH1*)eta_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][2] = (TH1*)deltaPhi_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][3] = (TH1*)helicityAngle_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][4] = (TH1*)pho1_pt_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][5] = (TH1*)pho2_pt_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][6] = (TH1*)pho1_eta_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][7] = (TH1*)pho2_eta_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][8] = (TH1*)pho_minr9_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][9] = (TH1*)maxeta_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][10] = (TH1*)ptOverM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][13] = (TH1*)ht_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][14] = (TH1*)deltaEta_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_GJetPF->Clone();
  hist_gjet_pf[3][21] = (TH1*)sigmaM_mhigh_cat2_GJetPF->Clone();

  hist2D_gjet_pf = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_GJetPF->Clone();
  hist_mass_gjet_pf = (TH1*)all_mass_cat0_GJetPF->Clone();
  hist_res_gjet_pf[0] = (TH1*)sigmaM_cat0_GJetPF->Clone();
  hist_res_gjet_pf[1] = (TH1*)sigmaM_cat1_GJetPF->Clone();
  hist_res_gjet_pf[2] = (TH1*)sigmaM_cat2_GJetPF->Clone();
  hist_res_gjet_pf[3] = (TH1*)sigmaM_cat3_GJetPF->Clone();
  hist_res_gjet_pf[4] = (TH1*)sigmaM_cat4_GJetPF->Clone();
  hist_res_gjet_pf[5] = (TH1*)sigmaM_cat5_GJetPF->Clone();
  hist_res_gjet_pf[6] = (TH1*)sigmaM_cat6_GJetPF->Clone();
  hist_res_gjet_pf[7] = (TH1*)sigmaM_cat7_GJetPF->Clone();
  hist_res_gjet_pf[8] = (TH1*)sigmaM_cat8_GJetPF->Clone();

  if (!madgraph) {

    f_qcd_pp->cd();
    hist_qcd_pp[1][0] = (TH1*)ptOverM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][1] = (TH1*)eta_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][2] = (TH1*)deltaPhi_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][3] = (TH1*)helicityAngle_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][4] = (TH1*)pho1_pt_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][5] = (TH1*)pho2_pt_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][6] = (TH1*)pho1_eta_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][7] = (TH1*)pho2_eta_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][8] = (TH1*)pho_minr9_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][9] = (TH1*)maxeta_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][10] = (TH1*)ptOverM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][13] = (TH1*)ht_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][14] = (TH1*)deltaEta_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][15] = (TH1*)deltaMOverMH_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][18] = (TH1*)sigmaMOverM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_QCDPP->Clone();
    hist_qcd_pp[1][21] = (TH1*)sigmaM_mlow_cat2_QCDPP->Clone();

    hist_qcd_pp[2][0] = (TH1*)ptOverM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][1] = (TH1*)eta_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][2] = (TH1*)deltaPhi_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][3] = (TH1*)helicityAngle_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][4] = (TH1*)pho1_pt_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][5] = (TH1*)pho2_pt_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][6] = (TH1*)pho1_eta_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][7] = (TH1*)pho2_eta_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][8] = (TH1*)pho_minr9_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][9] = (TH1*)maxeta_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][10] = (TH1*)ptOverM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][11] = (TH1*)pho1_ptOverM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][12] = (TH1*)pho2_ptOverM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][13] = (TH1*)ht_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][14] = (TH1*)deltaEta_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][15] = (TH1*)deltaMOverMH_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][18] = (TH1*)sigmaMOverM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_QCDPP->Clone();
    hist_qcd_pp[2][21] = (TH1*)sigmaM_msig_cat2_QCDPP->Clone();

    hist_qcd_pp[3][0] = (TH1*)ptOverM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][1] = (TH1*)eta_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][2] = (TH1*)deltaPhi_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][3] = (TH1*)helicityAngle_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][4] = (TH1*)pho1_pt_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][5] = (TH1*)pho2_pt_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][6] = (TH1*)pho1_eta_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][7] = (TH1*)pho2_eta_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][8] = (TH1*)pho_minr9_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][9] = (TH1*)maxeta_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][10] = (TH1*)ptOverM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][13] = (TH1*)ht_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][14] = (TH1*)deltaEta_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_QCDPP->Clone();
    hist_qcd_pp[3][21] = (TH1*)sigmaM_mhigh_cat2_QCDPP->Clone();

    hist2D_qcd_pp = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_QCDPP->Clone();
    hist_mass_qcd_pp = (TH1*)all_mass_cat0_QCDPP->Clone();
    hist_res_qcd_pp[0] = (TH1*)sigmaM_cat0_QCDPP->Clone();
    hist_res_qcd_pp[1] = (TH1*)sigmaM_cat1_QCDPP->Clone();
    hist_res_qcd_pp[2] = (TH1*)sigmaM_cat2_QCDPP->Clone();
    hist_res_qcd_pp[3] = (TH1*)sigmaM_cat3_QCDPP->Clone();
    hist_res_qcd_pp[4] = (TH1*)sigmaM_cat4_QCDPP->Clone();
    hist_res_qcd_pp[5] = (TH1*)sigmaM_cat5_QCDPP->Clone();
    hist_res_qcd_pp[6] = (TH1*)sigmaM_cat6_QCDPP->Clone();
    hist_res_qcd_pp[7] = (TH1*)sigmaM_cat7_QCDPP->Clone();
    hist_res_qcd_pp[8] = (TH1*)sigmaM_cat8_QCDPP->Clone();

  }

  f_qcd_pf->cd();
  hist_qcd_pf[1][0] = (TH1*)ptOverM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][1] = (TH1*)eta_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][2] = (TH1*)deltaPhi_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][3] = (TH1*)helicityAngle_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][4] = (TH1*)pho1_pt_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][5] = (TH1*)pho2_pt_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][6] = (TH1*)pho1_eta_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][7] = (TH1*)pho2_eta_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][8] = (TH1*)pho_minr9_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][9] = (TH1*)maxeta_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][10] = (TH1*)ptOverM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][13] = (TH1*)ht_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][14] = (TH1*)deltaEta_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][15] = (TH1*)deltaMOverMH_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][18] = (TH1*)sigmaMOverM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_QCDPF->Clone();
  hist_qcd_pf[1][21] = (TH1*)sigmaM_mlow_cat2_QCDPF->Clone();

  hist_qcd_pf[2][0] = (TH1*)ptOverM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][1] = (TH1*)eta_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][2] = (TH1*)deltaPhi_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][3] = (TH1*)helicityAngle_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][4] = (TH1*)pho1_pt_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][5] = (TH1*)pho2_pt_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][6] = (TH1*)pho1_eta_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][7] = (TH1*)pho2_eta_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][8] = (TH1*)pho_minr9_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][9] = (TH1*)maxeta_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][10] = (TH1*)ptOverM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][11] = (TH1*)pho1_ptOverM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][12] = (TH1*)pho2_ptOverM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][13] = (TH1*)ht_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][14] = (TH1*)deltaEta_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][15] = (TH1*)deltaMOverMH_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][18] = (TH1*)sigmaMOverM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_QCDPF->Clone();
  hist_qcd_pf[2][21] = (TH1*)sigmaM_msig_cat2_QCDPF->Clone();

  hist_qcd_pf[3][0] = (TH1*)ptOverM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][1] = (TH1*)eta_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][2] = (TH1*)deltaPhi_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][3] = (TH1*)helicityAngle_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][4] = (TH1*)pho1_pt_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][5] = (TH1*)pho2_pt_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][6] = (TH1*)pho1_eta_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][7] = (TH1*)pho2_eta_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][8] = (TH1*)pho_minr9_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][9] = (TH1*)maxeta_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][10] = (TH1*)ptOverM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][13] = (TH1*)ht_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][14] = (TH1*)deltaEta_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_QCDPF->Clone();
  hist_qcd_pf[3][21] = (TH1*)sigmaM_mhigh_cat2_QCDPF->Clone();

  hist2D_qcd_pf = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_QCDPF->Clone();
  hist_mass_qcd_pf = (TH1*)all_mass_cat0_QCDPF->Clone();
  hist_res_qcd_pf[0] = (TH1*)sigmaM_cat0_QCDPF->Clone();
  hist_res_qcd_pf[1] = (TH1*)sigmaM_cat1_QCDPF->Clone();
  hist_res_qcd_pf[2] = (TH1*)sigmaM_cat2_QCDPF->Clone();
  hist_res_qcd_pf[3] = (TH1*)sigmaM_cat3_QCDPF->Clone();
  hist_res_qcd_pf[4] = (TH1*)sigmaM_cat4_QCDPF->Clone();
  hist_res_qcd_pf[5] = (TH1*)sigmaM_cat5_QCDPF->Clone();
  hist_res_qcd_pf[6] = (TH1*)sigmaM_cat6_QCDPF->Clone();
  hist_res_qcd_pf[7] = (TH1*)sigmaM_cat7_QCDPF->Clone();
  hist_res_qcd_pf[8] = (TH1*)sigmaM_cat8_QCDPF->Clone();

  f_qcd_ff->cd();
  hist_qcd_ff[1][0] = (TH1*)ptOverM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][1] = (TH1*)eta_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][2] = (TH1*)deltaPhi_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][3] = (TH1*)helicityAngle_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][4] = (TH1*)pho1_pt_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][5] = (TH1*)pho2_pt_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][6] = (TH1*)pho1_eta_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][7] = (TH1*)pho2_eta_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][8] = (TH1*)pho_minr9_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][9] = (TH1*)maxeta_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][10] = (TH1*)ptOverM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][13] = (TH1*)ht_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][14] = (TH1*)deltaEta_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][15] = (TH1*)deltaMOverMH_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][18] = (TH1*)sigmaMOverM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_QCDFF->Clone();
  hist_qcd_ff[1][21] = (TH1*)sigmaM_mlow_cat2_QCDFF->Clone();

  hist_qcd_ff[2][0] = (TH1*)ptOverM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][1] = (TH1*)eta_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][2] = (TH1*)deltaPhi_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][3] = (TH1*)helicityAngle_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][4] = (TH1*)pho1_pt_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][5] = (TH1*)pho2_pt_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][6] = (TH1*)pho1_eta_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][7] = (TH1*)pho2_eta_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][8] = (TH1*)pho_minr9_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][9] = (TH1*)maxeta_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][10] = (TH1*)ptOverM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][11] = (TH1*)pho1_ptOverM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][12] = (TH1*)pho2_ptOverM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][13] = (TH1*)ht_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][14] = (TH1*)deltaEta_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][15] = (TH1*)deltaMOverMH_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][18] = (TH1*)sigmaMOverM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_QCDFF->Clone();
  hist_qcd_ff[2][21] = (TH1*)sigmaM_msig_cat2_QCDFF->Clone();

  hist_qcd_ff[3][0] = (TH1*)ptOverM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][1] = (TH1*)eta_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][2] = (TH1*)deltaPhi_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][3] = (TH1*)helicityAngle_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][4] = (TH1*)pho1_pt_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][5] = (TH1*)pho2_pt_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][6] = (TH1*)pho1_eta_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][7] = (TH1*)pho2_eta_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][8] = (TH1*)pho_minr9_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][9] = (TH1*)maxeta_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][10] = (TH1*)ptOverM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][13] = (TH1*)ht_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][14] = (TH1*)deltaEta_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_QCDFF->Clone();
  hist_qcd_ff[3][21] = (TH1*)sigmaM_mhigh_cat2_QCDFF->Clone();

  hist2D_qcd_ff = (TH2*)deltaMOverSigmaM_vs_sigmaM_cat2_QCDFF->Clone();
  hist_mass_qcd_ff = (TH1*)all_mass_cat0_QCDFF->Clone();
  hist_res_qcd_ff[0] = (TH1*)sigmaM_cat0_QCDFF->Clone();
  hist_res_qcd_ff[1] = (TH1*)sigmaM_cat1_QCDFF->Clone();
  hist_res_qcd_ff[2] = (TH1*)sigmaM_cat2_QCDFF->Clone();
  hist_res_qcd_ff[3] = (TH1*)sigmaM_cat3_QCDFF->Clone();
  hist_res_qcd_ff[4] = (TH1*)sigmaM_cat4_QCDFF->Clone();
  hist_res_qcd_ff[5] = (TH1*)sigmaM_cat5_QCDFF->Clone();
  hist_res_qcd_ff[6] = (TH1*)sigmaM_cat6_QCDFF->Clone();
  hist_res_qcd_ff[7] = (TH1*)sigmaM_cat7_QCDFF->Clone();
  hist_res_qcd_ff[8] = (TH1*)sigmaM_cat8_QCDFF->Clone();

  f_dy->cd();
  hist_dy[1][0] = (TH1*)ptOverM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][1] = (TH1*)eta_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][2] = (TH1*)deltaPhi_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][3] = (TH1*)helicityAngle_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][4] = (TH1*)pho1_pt_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][5] = (TH1*)pho2_pt_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][6] = (TH1*)pho1_eta_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][7] = (TH1*)pho2_eta_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][8] = (TH1*)pho_minr9_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][9] = (TH1*)maxeta_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][10] = (TH1*)ptOverM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][11] = (TH1*)pho1_ptOverM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][12] = (TH1*)pho2_ptOverM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][13] = (TH1*)ht_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][14] = (TH1*)deltaEta_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][15] = (TH1*)deltaMOverMH_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][16] = (TH1*)pho1_ptOverMH_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][17] = (TH1*)pho2_ptOverMH_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][18] = (TH1*)sigmaMOverM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][19] = (TH1*)deltaMSigmaMOverM2_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][20] = (TH1*)deltaMOverSigmaM_mlow_cat2_DyeeMG->Clone();
  hist_dy[1][21] = (TH1*)sigmaM_mlow_cat2_DyeeMG->Clone();

  hist_dy[2][0] = (TH1*)ptOverM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][1] = (TH1*)eta_msig_cat2_DyeeMG->Clone();
  hist_dy[2][2] = (TH1*)deltaPhi_msig_cat2_DyeeMG->Clone();
  hist_dy[2][3] = (TH1*)helicityAngle_msig_cat2_DyeeMG->Clone();
  hist_dy[2][4] = (TH1*)pho1_pt_msig_cat2_DyeeMG->Clone();
  hist_dy[2][5] = (TH1*)pho2_pt_msig_cat2_DyeeMG->Clone();
  hist_dy[2][6] = (TH1*)pho1_eta_msig_cat2_DyeeMG->Clone();
  hist_dy[2][7] = (TH1*)pho2_eta_msig_cat2_DyeeMG->Clone();
  hist_dy[2][8] = (TH1*)pho_minr9_msig_cat2_DyeeMG->Clone();
  hist_dy[2][9] = (TH1*)maxeta_msig_cat2_DyeeMG->Clone();
  hist_dy[2][10] = (TH1*)ptOverM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][11] = (TH1*)pho1_ptOverM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][12] = (TH1*)pho2_ptOverM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][13] = (TH1*)ht_msig_cat2_DyeeMG->Clone();
  hist_dy[2][14] = (TH1*)deltaEta_msig_cat2_DyeeMG->Clone();
  hist_dy[2][15] = (TH1*)deltaMOverMH_msig_cat2_DyeeMG->Clone();
  hist_dy[2][16] = (TH1*)pho1_ptOverMH_msig_cat2_DyeeMG->Clone();
  hist_dy[2][17] = (TH1*)pho2_ptOverMH_msig_cat2_DyeeMG->Clone();
  hist_dy[2][18] = (TH1*)sigmaMOverM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][19] = (TH1*)deltaMSigmaMOverM2_msig_cat2_DyeeMG->Clone();
  hist_dy[2][20] = (TH1*)deltaMOverSigmaM_msig_cat2_DyeeMG->Clone();
  hist_dy[2][21] = (TH1*)sigmaM_msig_cat2_DyeeMG->Clone();

  hist_dy[3][0] = (TH1*)ptOverM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][1] = (TH1*)eta_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][2] = (TH1*)deltaPhi_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][3] = (TH1*)helicityAngle_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][4] = (TH1*)pho1_pt_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][5] = (TH1*)pho2_pt_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][6] = (TH1*)pho1_eta_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][7] = (TH1*)pho2_eta_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][8] = (TH1*)pho_minr9_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][9] = (TH1*)maxeta_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][10] = (TH1*)ptOverM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][11] = (TH1*)pho1_ptOverM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][12] = (TH1*)pho2_ptOverM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][13] = (TH1*)ht_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][14] = (TH1*)deltaEta_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][15] = (TH1*)deltaMOverMH_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][16] = (TH1*)pho1_ptOverMH_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][17] = (TH1*)pho2_ptOverMH_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][18] = (TH1*)sigmaMOverM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][19] = (TH1*)deltaMSigmaMOverM2_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][20] = (TH1*)deltaMOverSigmaM_mhigh_cat2_DyeeMG->Clone();
  hist_dy[3][21] = (TH1*)sigmaM_mhigh_cat2_DyeeMG->Clone();

  hist_mass_dy = (TH1*)all_mass_cat0_DyeeMG->Clone();
  hist_res_dy[0] = (TH1*)sigmaM_cat0_DyeeMG->Clone();
  hist_res_dy[1] = (TH1*)sigmaM_cat1_DyeeMG->Clone();
  hist_res_dy[2] = (TH1*)sigmaM_cat2_DyeeMG->Clone();
  hist_res_dy[3] = (TH1*)sigmaM_cat3_DyeeMG->Clone();
  hist_res_dy[4] = (TH1*)sigmaM_cat4_DyeeMG->Clone();
  hist_res_dy[5] = (TH1*)sigmaM_cat5_DyeeMG->Clone();
  hist_res_dy[6] = (TH1*)sigmaM_cat6_DyeeMG->Clone();
  hist_res_dy[7] = (TH1*)sigmaM_cat7_DyeeMG->Clone();
  hist_res_dy[8] = (TH1*)sigmaM_cat8_DyeeMG->Clone();

  TString var[23] = {"ptOverM","eta","deltaPhi","helicityAngle","pho1_pt","pho2_pt","pho1_eta","pho2_eta","pho_minr9","maxeta","ptOverM","pho1_ptOverM","pho2_ptOverM","ht","deltaEta","deltaMOverMH","pho1_ptOverMH","pho2_ptOverMH","sigmaMOverM","deltaMSigmaMOverM2","deltaMOverSigmaM","sigmaM","bdtOut_grad"};
  TString title[23] = {"Diphoton p_{T}/M_{#gamma#gamma}","Diphoton #eta","#Delta#phi","cos#theta*","lead p_{T} (GeV)","sublead p_{T} (GeV)","lead #eta","sublead #eta","min(R9)","max #eta","Diphoton p_{T}/M_{#gamma#gamma}","p_{T1}/M_{#gamma#gamma}","p_{T2}/M_{#gamma#gamma}","H_T (GeV)","#Delta#eta","#DeltaM/M_{H}","p_{T1}/M_{H}","p_{T2}/M_{H}","#sigma_{M}/M_{#gamma#gamma}","#DeltaM/M_{H} * #sigma_{M}/M_{#gamma#gamma}","#DeltaM/#sigma_{M}","#sigma_{M}","BDT output (grad)"};
  if (equalBinWidths) {
    title[22]="BDT output bin number (gradient boost)";
    var[22]="bdtOutBin_grad";
    if (rebinBdtOut) var[22]="bdtOutBin2_grad";
  }

  TCanvas* canvas[23];

  float bornSF = madgraph ? 1.15 : 1.3;

  TCanvas *c_mgg2 = new TCanvas("c_mgg2","Mgg graph");
  c_mgg2->SetFillColor(0);

  float mass_hypothesis = float(mass_in);
  float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
  float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);

  //define sidebands
  float sideband_boundaries[4];
  sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
  sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
  sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
  sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);

  float Ndata_sig = hist_bkgModel[22]->Integral();
  float Ndata_sig_err;
  if (mass_in==115) { Ndata_sig_err=0.014; }
  else if (mass_in==120) { Ndata_sig_err=0.016; }
  else if (mass_in==121) { Ndata_sig_err=0.015; }
  else if (mass_in==123) { Ndata_sig_err=0.016; }
  else if (mass_in==125) { Ndata_sig_err=0.016; }
  else if (mass_in==130) { Ndata_sig_err=0.015; }
  else if (mass_in==135) { Ndata_sig_err=0.014; }
  else if (mass_in==140) { Ndata_sig_err=0.012; }
  else if (mass_in==150) { Ndata_sig_err=0.008; }
  else {
    cout << "Invalid mass: " << mass_in << endl;
  }

  for (int ivar=0; ivar<1; ivar++) {

    if (ivar==4 || ivar==5 || ivar==10 || ivar==11 || ivar==12 || ivar==13) continue;

    if (ivar==0 || ivar==10 || ivar==3) {
      hist_sig[ivar]->Rebin(2);
      for (int i=1; i<4; i++) {
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
      for (int i=1; i<4; i++) {
	hist_born[i][ivar]->Scale(1./0.6767);
	hist_box[i][ivar]->Scale(1./0.6767);
	if (!madgraph) {
	  hist_gjet_pp[i][ivar]->Scale(1./0.6767);
	  hist_qcd_pp[i][ivar]->Scale(1./0.6767);
	}
      }
    }

    if (data) hist_data[0][ivar] = (TH1*)hist_data[1][ivar]->Clone();
    hist_born[0][ivar] = (TH1*)hist_born[1][ivar]->Clone();
    hist_box[0][ivar] = (TH1*)hist_box[1][ivar]->Clone();
    hist_gjet_pf[0][ivar] = (TH1*)hist_gjet_pf[1][ivar]->Clone();
    hist_qcd_pf[0][ivar] = (TH1*)hist_qcd_pf[1][ivar]->Clone();
    hist_qcd_ff[0][ivar] = (TH1*)hist_qcd_ff[1][ivar]->Clone();
    hist_dy[0][ivar] = (TH1*)hist_dy[1][ivar]->Clone();
    if (!madgraph) {
      hist_gjet_pp[0][ivar] = (TH1*)hist_gjet_pp[1][ivar]->Clone();
      hist_qcd_pp[0][ivar] = (TH1*)hist_qcd_pp[1][ivar]->Clone();
    }
    for (int i=2; i<4; i++) {
      if (data) hist_data[0][ivar]->Add(hist_data[i][ivar]);
      hist_born[0][ivar]->Add(hist_born[i][ivar]);
      hist_box[0][ivar]->Add(hist_box[i][ivar]);
      hist_gjet_pf[0][ivar]->Add(hist_gjet_pf[i][ivar]);
      hist_qcd_pf[0][ivar]->Add(hist_qcd_pf[i][ivar]);
      hist_qcd_ff[0][ivar]->Add(hist_qcd_ff[i][ivar]);
      hist_dy[0][ivar]->Add(hist_dy[i][ivar]);
      if (!madgraph) {
	hist_gjet_pp[0][ivar]->Add(hist_gjet_pp[i][ivar]);
	hist_qcd_pp[0][ivar]->Add(hist_qcd_pp[i][ivar]);
      }
    }

    //Apply k-factors
    for (int i=0; i<4; i++) {
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

    //Apply +ve/-ve folding for eta distributions
    if (ivar==1 || ivar==6 ||ivar==7 || ivar==9) {
      for (int ibin=25; ibin<50; ibin++) {
	hist_sig[ivar]->SetBinContent(ibin+1,hist_sig[ivar]->GetBinContent(ibin+1)+hist_sig[ivar]->GetBinContent(50-ibin));
	for (int i=0; i<4; i++) {
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
	for (int i=0; i<4; i++) {
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

    for (int i=0; i<4; i++) {
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
    if (data) hist_data[0][ivar]->SetMarkerStyle(20);
    if (data) hist_data[0][ivar]->SetMarkerSize(.8);
    if (data) hist_data[2][ivar]->SetMarkerStyle(20);
    if (data) hist_data[2][ivar]->SetMarkerSize(.8);
    for (int i=0; i<4; i++) {
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
    hist_bkg_stack[ivar]->Add(hist_box[0][ivar]);
    hist_bkg_stack[ivar]->Add(hist_born[0][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_box[2][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_born[2][ivar]);
    if (!madgraph) {
      hist_bkg_stack[ivar]->Add(hist_gjet_pp[0][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_pp[0][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_gjet_pp[2][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_pp[2][ivar]);
    }
    if (fakes) {
      hist_bkg_stack[ivar]->Add(hist_gjet_pf[0][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_pf[0][ivar]);
      hist_bkg_stack[ivar]->Add(hist_qcd_ff[0][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_gjet_pf[2][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_pf[2][ivar]);
      hist_bkg_stack_sig[ivar]->Add(hist_qcd_ff[2][ivar]);
    }    
    hist_bkg_stack[ivar]->Add(hist_dy[0][ivar]);
    hist_bkg_stack_sig[ivar]->Add(hist_dy[2][ivar]);

    for (int i=0; i<4; i++) {
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
    hist_sig[ivar]->Scale(hist_bkg[2][ivar]->Integral()/nsig);

    canvas[ivar] = new TCanvas("c_"+var[ivar],var[ivar],1950,1300);
    canvas[ivar]->Divide(3,3);
    canvas[ivar]->SetFillColor(0);

    canvas[ivar]->cd(1);

    TLegend *leg;
    if (ivar==0 || ivar==1 || ivar==4 || ivar==5 || ivar==6 || ivar==7 || (ivar>=10 && ivar!=14)) {
      leg = new TLegend(.55,.55,.87,.87);
    } else if (ivar==3 || ivar==14) {
      leg = new TLegend(.32,.1,.68,.35);
    } else {
      leg = new TLegend(.15,.55,.47,.87);
    }
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetTextSize(.035);
    if (data) leg->AddEntry(hist_data[2][ivar],"Data (1.66fb^{-1})");
    leg->AddEntry(hist_sig[ivar],"GluGlu"+mass_str);
    leg->AddEntry(hist_dy[2][ivar],"DYee+Z","F");
    if (fakes) {
      leg->AddEntry(hist_qcd_ff[2][ivar],"QCD fake-fake","F");
      leg->AddEntry(hist_qcd_pf[2][ivar],"QCD prompt-fake","F");
      leg->AddEntry(hist_gjet_pf[2][ivar],"GJet prompt-fake","F");
    }
    if (!madgraph) {
      leg->AddEntry(hist_qcd_pp[2][ivar],"QCD prompt-prompt","F");
      leg->AddEntry(hist_gjet_pp[2][ivar],"GJet prompt-prompt","F");
    }
    if (!madgraph) {
      leg->AddEntry(hist_born[2][ivar],"Pythia Born","F");
    } else {
      leg->AddEntry(hist_born[2][ivar],"Madgraph DiPhotonJets","F");
    }
    leg->AddEntry(hist_box[2][ivar],"Pythia Box","F");


    float max = hist_bkg_stack_sig[ivar]->GetMaximum();
    if (hist_sig[ivar]->GetMaximum()>max) max=hist_sig[ivar]->GetMaximum();
    if (data && hist_data[2][ivar]->GetMaximum()>max) max=hist_data[2][ivar]->GetMaximum();
    hist_bkg_stack_sig[ivar]->SetMaximum(max*1.1);

    hist_bkg_stack_sig[ivar]->Draw("hist");
    if (ivar==1) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    //if (ivar==21) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitleSize(0.04);
    if (ivar!=9 && ivar!=3 && ivar!=14) leg->Draw();
    hist_sig[ivar]->Draw("same");
    if (data) hist_data[2][ivar]->Draw("same,e");

    txt = new TLatex();
    txt->SetNDC();
    txt->SetTextSize(0.045);
    if (ivar==2 || ivar==8) {
      txt->DrawLatex(0.15,0.42, "Data and MC in signal region");
    } else {
      txt->DrawLatex(0.15,0.82, "Data and MC in signal region");
    }

    canvas[ivar]->cd(3);

    hist_bkg_reweight[ivar] = (TH1*)hist_bkg[1][ivar]->Clone();
    float nbkg_sig = hist_bkg[2][ivar]->Integral();
    hist_bkg_reweight[ivar]->Reset();
    hist_bkg_reweight[ivar]->Add(hist_bkg[1][ivar],0.5*nbkg_sig/hist_bkg[1][ivar]->Integral());
    hist_bkg_reweight[ivar]->Add(hist_bkg[3][ivar],0.5*nbkg_sig/hist_bkg[3][ivar]->Integral());

    int rebinFac = ivar==15 ? 5 : 2;
    if (ivar==1 || ivar==6 ||ivar==7 || ivar==9 || ivar==3) rebinFac=1;
    if (ivar==19) rebinFac=4;
    hist_bkg_reweight[ivar]->Rebin(rebinFac);
    float sf=hist_bkg_reweight[ivar]->Integral();
    hist_bkg_reweight[ivar]->Scale(1./sf);
    for (int i=1; i<4; i++) {
      hist_bkg[i][ivar]->Rebin(rebinFac);
      float sf=hist_bkg[i][ivar]->Integral();
      hist_bkg[i][ivar]->Scale(1./sf);
    }

    max = hist_bkg[1][ivar]->GetMaximum();
    if (hist_bkg[2][ivar]->GetMaximum()>max) max=hist_bkg[2][ivar]->GetMaximum();
    if (hist_bkg[3][ivar]->GetMaximum()>max) max=hist_bkg[3][ivar]->GetMaximum();
    hist_bkg[1][ivar]->SetMaximum(max*1.1);

    hist_bkg[1][ivar]->SetLineColor(2);
    hist_bkg[2][ivar]->SetLineColor(4);
    hist_bkg[3][ivar]->SetLineColor(6);
    hist_bkg[2][ivar]->SetLineWidth(2);
    hist_bkg[1][ivar]->SetLineWidth(2);
    hist_bkg[3][ivar]->SetLineWidth(2);
    //hist_bkg[1][ivar]->SetLineStyle(2);
    //hist_bkg[3][ivar]->SetLineStyle(3);
    hist_bkg[1][ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg[1][ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_bkg[1][ivar]->GetXaxis()->SetRangeUser(0.,7.);
    //if (ivar==21) hist_bkg[1][ivar]->GetXaxis()->SetRangeUser(0.,10.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg[1][ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    hist_bkg[1][ivar]->SetFillColor(2);
    //hist_bkg[1][ivar]->SetFillStyle(3004);
    hist_bkg[3][ivar]->SetFillColor(6);
    //hist_bkg[3][ivar]->SetFillStyle(3005);
    hist_bkg[2][ivar]->SetFillColor(4);
    hist_bkg[1][ivar]->Draw("e3");
    hist_bkg[3][ivar]->Draw("e3,same");
    hist_bkg[2][ivar]->Draw("e3,same");

    TLegend *leg2 = (TLegend*)leg->Clone();
    leg2->Clear();
    leg2->AddEntry(hist_bkg[3][ivar],"MC High sideband");
    leg2->AddEntry(hist_bkg[2][ivar],"MC Signal region");
    leg2->AddEntry(hist_bkg[1][ivar],"MC Low sideband");
    if (ivar!=9) leg2->Draw();

    canvas[ivar]->cd(2);

    hist_data_reweight[ivar] = (TH1*)hist_data[1][ivar]->Clone();
    hist_data_reweight[ivar]->Reset();
    hist_data_reweight[ivar]->Add(hist_data[1][ivar],0.5*Ndata_sig/hist_data[1][ivar]->Integral());
    hist_data_reweight[ivar]->Add(hist_data[3][ivar],0.5*Ndata_sig/hist_data[3][ivar]->Integral());
    //minr9 histogram is filled twice (once for each photon) so needs twice the normalization
    //if (ivar==8) hist_data_reweight[ivar]->Scale(2.);

    hist_sig_reweight[ivar]->Scale(nbkg_sig/nsig);

    max = hist_bkg_stack_sig[ivar]->GetMaximum();
    if (hist_sig_reweight[ivar]->GetMaximum()>max) max=hist_sig_reweight[ivar]->GetMaximum();
    if (data && hist_data_reweight[ivar]->GetMaximum()>max) max=hist_data_reweight[ivar]->GetMaximum();
    hist_bkg_stack_sig[ivar]->SetMaximum(max*1.1);

    hist_sig_reweight[ivar]->SetLineColor(4);
    hist_sig_reweight[ivar]->SetLineWidth(2.5);
    if (data) hist_data_reweight[ivar]->SetMarkerStyle(20);
    if (data) hist_data_reweight[ivar]->SetMarkerSize(.8);

    hist_bkg_stack_sig[ivar]->Draw("hist");
    if (ivar==1) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    //if (ivar==21) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg_stack_sig[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg_stack_sig[ivar]->GetXaxis()->SetTitleSize(0.04);
    if (ivar!=9 && ivar!=3 && ivar!=14) leg->Draw();
    hist_sig_reweight[ivar]->Draw("same");
    if (data) hist_data_reweight[ivar]->Draw("same,e");

    if (ivar==2 || ivar==8) {
      txt->DrawLatex(0.15,0.42,"Data from weighted sidebands");
      txt->DrawLatex(0.15,0.36,"MC in signal region");
    } else {
      txt->DrawLatex(0.15,0.82,"Data from weighted sidebands");
      txt->DrawLatex(0.15,0.76,"MC in signal region");
    }

    canvas[ivar]->cd(4);

    hist_bkg_sig[ivar] = (TH1*)hist_bkg[2][ivar]->Clone();

    max = hist_bkg_reweight[ivar]->GetMaximum();
    if (hist_bkg_sig[ivar]->GetMaximum()>max) max=hist_bkg_sig[ivar]->GetMaximum();
    hist_bkg_reweight[ivar]->SetMaximum(max*1.1);

    hist_bkg_reweight[ivar]->SetMarkerStyle(20);
    hist_bkg_reweight[ivar]->SetMarkerSize(.8);
    hist_bkg_reweight[ivar]->SetMarkerColor(4);
    hist_bkg_reweight[ivar]->SetLineColor(4);
    hist_bkg_reweight[ivar]->SetLineWidth(2);
    hist_bkg_reweight[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_bkg_reweight[ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    //if (ivar==21) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_bkg_reweight[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    hist_bkg_sig[ivar]->SetFillColor(38);
    hist_bkg_reweight[ivar]->Draw("e");
    hist_bkg_sig[ivar]->Draw("e2,same");
    hist_bkg_reweight[ivar]->Draw("e,same");

    TLegend *leg3 = (TLegend*)leg->Clone();
    leg3->Clear();
    leg3->AddEntry(hist_bkg_sig[ivar],"MC Signal region");
    leg3->AddEntry(hist_bkg_reweight[ivar],"MC Background model");
    if (ivar!=9) leg3->Draw();

    canvas[ivar]->cd(5);

    hist_data_sig[ivar] = (TH1*)hist_data[2][ivar]->Clone();
    hist_data_sig[ivar]->SetMarkerStyle(0);
    hist_data_sig[ivar]->Scale(hist_data_reweight[ivar]->Integral()/hist_data_sig[ivar]->Integral());

    hist_data_reweight_rebin[ivar] = (TH1*)hist_data_reweight[ivar]->Clone();
    hist_data_reweight_rebin[ivar]->Rebin(rebinFac);
    hist_data_sig[ivar]->Rebin(rebinFac);

    max = hist_data_reweight_rebin[ivar]->GetMaximum();
    if (hist_data_sig[ivar]->GetMaximum()>max) max=hist_data_sig[ivar]->GetMaximum();
    hist_data_reweight_rebin[ivar]->SetMaximum(max*1.1);

    hist_data_reweight_rebin[ivar]->GetXaxis()->SetTitle(title[ivar]);
    hist_data_reweight_rebin[ivar]->GetXaxis()->SetTitleSize(0.04);

    if (ivar==1) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,7.);
    //if (ivar==21) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,10.);
    if (ivar==6 || ivar==7 ||ivar==9) hist_data_reweight_rebin[ivar]->GetXaxis()->SetRangeUser(0.,2.5);
    hist_data_sig[ivar]->SetFillColor(38);
    hist_data_reweight_rebin[ivar]->Draw("e");
    hist_data_sig[ivar]->Draw("e2,same");
    hist_data_reweight_rebin[ivar]->Draw("e,same");

    TLegend *leg3 = (TLegend*)leg->Clone();
    leg3->Clear();
    leg3->AddEntry(hist_data_sig[ivar],"Data: Signal region");
    leg3->AddEntry(hist_data_reweight_rebin[ivar],"Background Model");
    if (ivar!=9) leg3->Draw();

    canvas[ivar]->cd(7);

    hist_bkg_ratio[ivar] = (TH1*)hist_bkg_reweight[ivar]->Clone();
    hist_bkg_ratio[ivar]->Reset();
    hist_bkg_ratio[ivar]->Add(hist_bkg_reweight[ivar]);
    hist_bkg_ratio[ivar]->Divide(hist_bkg_sig[ivar]);

    hist_bkg_ratio[ivar]->SetMarkerStyle(20);
    hist_bkg_ratio[ivar]->SetMarkerColor(4);
    hist_bkg_ratio[ivar]->SetMarkerSize(.8);

    hist_bkg_ratio[ivar]->SetMaximum(1.4);
    hist_bkg_ratio[ivar]->SetMinimum(0.6);
    hist_bkg_ratio[ivar]->Draw();

    float xmin = hist_bkg_ratio[ivar]->GetXaxis()->GetXmin();
    if (ivar==1 || ivar==6 || ivar==7 ||ivar==9) xmin=0.;
    float xmax = hist_bkg_ratio[ivar]->GetXaxis()->GetXmax();
    TLine *line2 = new TLine(xmin,1.,xmax,1.);
    line2->SetLineColor(4);
    line2->SetLineWidth(2);
    line2->Draw();

    txt->DrawLatex(0.15,0.82, "MC weighted sidebands / MC in signal region");

    canvas[ivar]->cd(8);
    
    hist_data_ratio[ivar] = (TH1*)hist_data_reweight_rebin[ivar]->Clone();
    hist_data_ratio[ivar]->Reset();
    hist_data_ratio[ivar]->Add(hist_data_reweight_rebin[ivar]);
    hist_data_ratio[ivar]->Divide(hist_data_sig[ivar]);

    hist_data_ratio[ivar]->SetMarkerStyle(20);
    hist_data_ratio[ivar]->SetMarkerSize(.8);

    hist_data_ratio[ivar]->SetMaximum(1.4);
    hist_data_ratio[ivar]->SetMinimum(0.6);
    hist_data_ratio[ivar]->Draw("e");

    line2->Draw();

    txt->DrawLatex(0.15,0.82,"Background Model / Data in signal region");

    canvas[ivar]->SaveAs(outdir+var[ivar]+".gif");

  }

  if (mass_in != 110) {

    hist_sig[22]->SetLineColor(2);
    hist_sig[22]->SetLineWidth(2.5);

    for (int i=1; i<4; i++) {
      hist_data[i][22]->SetMarkerStyle(20);
      hist_data[i][22]->SetMarkerSize(.8);
      hist_data[i][22]->SetLineWidth(2);
      hist_data[i][22]->GetXaxis()->SetTitle(title[22]);
      hist_data[i][22]->GetXaxis()->SetTitleSize(0.04);
      hist_bkg[i][22]->SetLineWidth(2);
      hist_bkg[i][22]->GetXaxis()->SetTitle(title[22]);
      hist_bkg[i][22]->GetXaxis()->SetTitleSize(0.04);
    }

    hist_bkg[1][22]->SetLineColor(2);
    hist_bkg[1][22]->SetMarkerColor(2);
    hist_bkg[1][22]->SetMarkerStyle(20);
    hist_bkg[1][22]->SetMarkerSize(0.8);
    hist_bkg[2][22]->SetLineColor(4);
    hist_bkg[2][22]->SetMarkerColor(4);
    hist_bkg[2][22]->SetFillColor(38);
    hist_bkgModel[22]->SetMarkerStyle(20);
    hist_bkgModel[22]->SetMarkerSize(.8);
    hist_bkgModel[22]->SetLineColor(4);
    hist_bkgModel[22]->SetFillColor(38);
    hist_bkgModel[22]->SetLineWidth(2);
    hist_bkgModel[22]->GetXaxis()->SetTitle(title[22]);
    hist_bkgModel[22]->GetXaxis()->SetTitleSize(0.04);
    hist_bkg[3][22]->SetLineColor(kGreen-2);
    hist_bkg[3][22]->SetMarkerColor(kGreen-2);
    hist_bkg[3][22]->SetMarkerStyle(20);
    hist_bkg[3][22]->SetMarkerSize(0.8);
    hist_sig_x5 = (TH1*)hist_sig[22]->Clone();
    hist_sig_x10 = (TH1*)hist_sig[22]->Clone();
    hist_sig_x5->SetLineStyle(2);
    hist_sig[22]->SetLineStyle(3);
    hist_sig[22]->Scale(2.);
    hist_sig_x5->Scale(5.);
    hist_sig_x10->Scale(10.);

    for (int i=1; i<4; i++) {
      float sf=nbkg_sig/hist_bkg[i][22]->Integral();
      hist_bkg_scaled[i][22] = (TH1*)hist_bkg[i][22]->Clone();
      hist_bkg_scaled[i][22]->Scale(sf);

      sf=Ndata_sig/hist_data[i][22]->Integral();
      hist_data_scaled[i][22] = (TH1*)hist_data[i][22]->Clone();
      hist_data_scaled[i][22]->Scale(sf);

    }

    hist_bkg_reweight[22] = (TH1*)hist_bkg[1][22]->Clone();
    float nbkg_sig = hist_bkg[2][22]->Integral();
    hist_bkg_reweight[22]->Reset();
    hist_bkg_reweight[22]->Add(hist_bkg[1][22],0.5*nbkg_sig/hist_bkg[1][22]->Integral());
    hist_bkg_reweight[22]->Add(hist_bkg[3][22],0.5*nbkg_sig/hist_bkg[3][22]->Integral());

    canvas[22] = new TCanvas("c_"+var[22],var[22],1950,1300);
    canvas[22]->Divide(3,3);
    canvas[22]->SetFillColor(0);

    canvas[22]->cd(1);

    hist_bkgModel[22]->SetMaximum(150.);
    hist_bkgModel[22]->SetMinimum(0.);

    hist_bkgModel[22]->Draw("hist");
    /*
    hist_bkgModel_err = (TH1*)hist_bkgModel[22]->Clone();
    for (int ibin=0; ibin<nbins+1; ibin++) {
      //hist_bkgModel_err->SetBinError(ibin,hist_bkgModel[22]->GetBinContent(ibin)*Ndata_sig_err);
    }
    hist_bkgModel_err->SetFillColor(4);
    hist_bkgModel_err->Draw("e2,same");
    */
    hist_sig[22]->Draw("hist,same");
    hist_sig_x5->Draw("hist,same");
    hist_sig_x10->Draw("hist,same");
    hist_data[2][22]->Draw("same,e");

    TLegend *leg1 = new TLegend(.15,.55,.47,.87);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(10);
    leg1->SetTextSize(.035);
    leg1->AddEntry(hist_data[2][22],"Data (1.66fb^{-1})");
    leg1->AddEntry(hist_sig_x10,"Signal ("+mass_str+" GeV) x2, x5, x10");
    leg1->AddEntry(hist_bkgModel[22],"Background Model","F");
    leg1->Draw();
    leg1->Draw();

    canvas[22]->cd(5);

    hist_data_highMinusLow = (TH1*)hist_data_scaled[3][22]->Clone();
    hist_data_highMinusLow->Add(hist_data_scaled[1][22],-1.);
    hist_data_highMinusLow->SetMarkerColor(4);
    hist_data_highMinusLow->SetLineColor(4);
    hist_data_highMinusLow->SetMaximum(50.);
    hist_data_highMinusLow->SetMinimum(-50.);
    hist_data_highMinusLow->Draw("e");

    txt->DrawLatex(0.15,0.82, "high sideband - low sideband");

    float xmin = hist_data[1][22]->GetXaxis()->GetXmin();
    float xmax = hist_data[1][22]->GetXaxis()->GetXmax();
    TLine *line1 = new TLine(xmin,0.,xmax,0.);
    line1->SetLineColor(4);
    line1->SetLineWidth(2);
    line1->Draw();

    /*
    canvas[22]->cd(3);

    hist_bkg[2][22]->SetMaximum(150.);
    hist_bkg[2][22]->SetMinimum(0.);

    hist_bkg[2][22]->Draw("hist");
    hist_sig[22]->Draw("hist,same");
    hist_data[2][22]->Draw("same,e");

    leg3 = (TLegend*)leg1->Clone();
    leg3->Clear();
    leg3->AddEntry(hist_data[2][22],"Data (1.66fb^{-1})");
    leg3->AddEntry(hist_sig[22],"Signal ("+mass_str+" GeV) x 10");
    leg3->AddEntry(hist_bkg[2][22],"Background MC","F");
    leg3->Draw();
    */

    canvas[22]->cd(6);

    hist_bkg_highMinusLow = (TH1*)hist_bkg_scaled[3][22]->Clone();
    hist_bkg_highMinusLow->Add(hist_bkg_scaled[1][22],-1.);
    hist_bkg_highMinusLow->SetMarkerColor(4);
    hist_bkg_highMinusLow->SetLineColor(4);
    hist_bkg_highMinusLow->SetMaximum(50.);
    hist_bkg_highMinusLow->SetMinimum(-50.);
    hist_bkg_highMinusLow->Draw("e");

    txt->DrawLatex(0.15,0.82, "MC high sideband - MC low sideband");

    line1->Draw();

    canvas[22]->cd(4);

    hist_dataMinusModel = (TH1*)hist_data[2][22]->Clone();
    hist_dataMinusModel->Add(hist_bkgModel[22],-1.);

    float max=50.;
    if (hist_sig[22]->GetMaximum()>max) max=hist_sig[22]->GetMaximum()*1.05;
    hist_dataMinusModel->SetMaximum(max);
    hist_dataMinusModel->SetMinimum(-50.);
    hist_dataMinusModel->Draw("e");
    hist_sig[22]->Draw("hist,same");
    hist_sig_x5->Draw("hist,same");
    hist_sig_x10->Draw("hist,same");

    TLegend *leg2 = new TLegend(.15,.7,.47,.87);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(10);
    leg2->SetTextSize(.035);
    leg2->AddEntry(hist_dataMinusModel,"Data (1.66fb^{-1}) - background model","LP");
    leg2->AddEntry(hist_sig_x10,"Signal ("+mass_str+" GeV) x2, x5, x10");
    leg2->Draw();

    line1->Draw();

    canvas[22]->cd(2);

    hist_data_scaled[1][22]->SetLineColor(2);
    hist_data_scaled[1][22]->SetMarkerColor(2);
    hist_data_scaled[3][22]->SetLineColor(kGreen-2);
    hist_data_scaled[3][22]->SetMarkerColor(kGreen-2);

    hist_data_scaled[1][22]->SetMaximum(150.);
    hist_data_scaled[1][22]->SetMinimum(0.);
    hist_data_scaled[1][22]->Draw("e");
    hist_data_scaled[3][22]->Draw("e,same");

    leg6 = (TLegend*)leg2->Clone();
    leg6->Clear();
    leg6->AddEntry(hist_data_scaled[3][22],"High sideband","LP");
    leg6->AddEntry(hist_data_scaled[1][22],"Low sideband","LP");
    leg6->Draw();

    canvas[22]->cd(9);

    hist_bkg_reweight[22]->SetMaximum(150.);
    hist_bkg_reweight[22]->SetMinimum(0.);
    hist_bkg_reweight[22]->SetLineColor(4);
    hist_bkg_reweight[22]->SetMarkerColor(4);
    hist_bkg[2][22]->SetFillColor(38);

    hist_bkg_reweight[22]->Draw("e");
    hist_bkg[2][22]->Draw("e2,same");
    hist_bkg_reweight[22]->Draw("e,same");

    TLegend *leg7 = (TLegend*)leg1->Clone();
    leg7->Clear();
    leg7->AddEntry(hist_bkg[2][22],"MC Signal region","F");
    leg7->AddEntry(hist_bkg_reweight[22],"MC Background model");
    leg7->Draw();

    canvas[22]->cd(3);

    hist_bkg_scaled[2][22]->SetMaximum(150.);
    hist_bkg_scaled[2][22]->SetMinimum(0.);

    hist_bkg_scaled[2][22]->Draw("e2");
    hist_bkg_scaled[1][22]->Draw("e,same");
    hist_bkg_scaled[3][22]->Draw("e,same");

    leg8 = (TLegend*)leg1->Clone();
    leg8->Clear();
    leg8->AddEntry(hist_bkg_scaled[3][22],"MC High sideband","LP");
    leg8->AddEntry(hist_bkg_scaled[2][22],"MC Signal region","F");
    leg8->AddEntry(hist_bkg_scaled[1][22],"MC Low sideband","LP");
    leg8->Draw();

    /*
    canvas[22]->cd(11);

    hist_modelMinusSig = (TH1*)hist_bkg_reweight[22]->Clone();
    hist_modelMinusSig->Add(hist_bkg[2][22],-1.);

    hist_modelMinusSig->SetMaximum(50.);
    hist_modelMinusSig->SetMinimum(-50.);
    hist_modelMinusSig->Draw();

    line1->Draw();

    txt->DrawLatex(0.15,0.82, "MC Background model - MC signal region");

    canvas[22]->cd(12);

    hist_bkg_highMinusSig = (TH1*)hist_bkg_scaled[3][22]->Clone();
    hist_bkg_lowMinusSig = (TH1*)hist_bkg_scaled[1][22]->Clone();
    hist_bkg_highMinusSig->Add(hist_bkg_scaled[2][22],-1.);
    hist_bkg_lowMinusSig->Add(hist_bkg_scaled[2][22],-1.);

    hist_bkg_highMinusSig->SetMarkerColor(kGreen-2);
    hist_bkg_highMinusSig->SetLineColor(kGreen-2);
    hist_bkg_lowMinusSig->SetMarkerColor(2);
    hist_bkg_lowMinusSig->SetLineColor(2);

    hist_bkg_highMinusSig->SetMaximum(50.);
    hist_bkg_highMinusSig->SetMinimum(-50.);
    hist_bkg_highMinusSig->Draw("e");
    hist_bkg_lowMinusSig->Draw("e,same");

    leg12 = (TLegend*)leg2->Clone();
    leg12->Clear();
    leg12->AddEntry(hist_bkg_highMinusSig,"MC high sideband - MC signal region");
    leg12->AddEntry(hist_bkg_lowMinusSig,"MC high sideband - MC signal region");
    leg12->Draw();

    line1->Draw();
    */

    canvas[22]->SaveAs(outdir+var[22]+".gif");

  }

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

  cout << "Data/MC scale factor = " << hist_mass_data->Integral()/hist_mass_bkg->Integral() << endl;

  TCanvas *c_mgg = new TCanvas("c_mgg","Mgg",1000,700);
  c_mgg->SetFillColor(0);

  float max = hist_mass_bkg_stack->GetMaximum();
  if (hist_mass_sig->GetMaximum()>max) max=hist_mass_sig->GetMaximum();
  if (data && hist_mass_data->GetMaximum()>max) max=hist_mass_data->GetMaximum();
  //max=380.;
  hist_mass_bkg_stack->SetMaximum(max*1.05);

  hist_mass_bkg_stack->Draw();
  //hist_mass_bkg_stack->GetXaxis()->SetRangeUser(95.,185.);
  hist_mass_bkg_stack->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  hist_mass_bkg_stack->GetXaxis()->SetTitleSize(0.04);
  hist_mass_sig->Draw("same");
  if (data) hist_mass_data->Draw("same,e");

  TLine* line[4];
  for (int i=0; i<4; i++) {
    line[i] = new TLine(sideband_boundaries[i],0.,sideband_boundaries[i],max*1.1);
    line[i]->SetLineColor(38);
    line[i]->SetLineWidth(3);
    line[i]->SetLineStyle(9);
    line[i]->Draw();
  }

  TLegend *leg_mass = new TLegend(.65,.6,.87,.87);
  leg_mass->SetBorderSize(0);
  leg_mass->SetFillColor(10);
  leg_mass->SetTextSize(.025);
  if (data) leg_mass->AddEntry(hist_mass_data,"Data (1.66fb^{-1})");
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

    hist_res_born[icat]->Scale(bornSF*dataMC_sf);
    hist_res_box[icat]->Scale(1.3*dataMC_sf);
    hist_res_gjet_pf[icat]->Scale(1.3*dataMC_sf);
    hist_res_qcd_pf[icat]->Scale(1.3*dataMC_sf);
    hist_res_qcd_ff[icat]->Scale(1.*dataMC_sf);
    hist_res_dy[icat]->Scale(1.15*2321./992.*dataMC_sf);
    if (!madgraph) {
      hist_res_gjet_pp[icat]->Scale(1.3*dataMC_sf);
      hist_res_qcd_pp[icat]->Scale(1.3*dataMC_sf);
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

  //TString title2D[2] = {"M_{#gamma#gamma}","H_{T}"};
  //TString title2D[2] = {"#DeltaM/M_{H}","#sigma_{M}"};
  //TString title2D[2] = {"#DeltaM/M_{H}","#DeltaM/Sigma_{M}"};
  TString title2D[2] = {"#sigma_{M}","#DeltaM/Sigma_{M}"};

  /*
  hist2D_born->Scale(bornSF);
  hist2D_box->Scale(1.3);
  hist2D_gjet_pf->Scale(1.3);
  hist2D_qcd_pf->Scale(1.3);
  hist2D_qcd_ff->Scale(1.);
  if (!madgraph) {
    hist2D_gjet_pp->Scale(1.3);
    hist2D_qcd_pp->Scale(1.3);
  }

  hist2D_bkg = (TH2*)hist2D_box->Clone();
  hist2D_bkg->Add(hist2D_born);
  if (!madgraph) {
    hist2D_bkg->Add(hist2D_gjet_pp);
    hist2D_bkg->Add(hist2D_qcd_pp);
  }
  if (fakes) {
    hist2D_bkg->Add(hist2D_gjet_pf);
    hist2D_bkg->Add(hist2D_qcd_pf);
    hist2D_bkg->Add(hist2D_qcd_ff);
  }
  hist2D_bkg->GetXaxis()->SetTitle(title2D[0]);
  hist2D_bkg->GetYaxis()->SetTitle(title2D[1]);
  hist2D_bkg->GetXaxis()->SetTitleSize(0.04);
  hist2D_bkg->GetYaxis()->SetTitleSize(0.04);
  hist2D_bkg->SetFillColor(kGreen-2);

  hist2D_sig->GetXaxis()->SetTitle(title2D[0]);
  hist2D_sig->GetYaxis()->SetTitle(title2D[1]);
  hist2D_sig->GetXaxis()->SetTitleSize(0.04);
  hist2D_sig->GetYaxis()->SetTitleSize(0.04);
  hist2D_sig->SetFillColor(kBlue);

  if (data) {
    hist2D_data->GetXaxis()->SetTitle(title2D[0]);
    hist2D_data->GetYaxis()->SetTitle(title2D[1]);
    hist2D_data->GetXaxis()->SetTitleSize(0.04);
    hist2D_data->GetYaxis()->SetTitleSize(0.04);
  }

  TCanvas *c_deltaMOverSigmaM_vs_sigmaM = new TCanvas("c_deltaMOverSigmaM_vs_sigmaM","deltaMOverSigmaM_vs_sigmaM",1200,800);
  c_deltaMOverSigmaM_vs_sigmaM->Divide(2,2);
  c_deltaMOverSigmaM_vs_sigmaM->SetFillColor(10);
  c_deltaMOverSigmaM_vs_sigmaM->cd(1);
  hist2D_bkg->Draw("box");
  txt->DrawLatex(0.15,0.82,"Total Background");
  c_deltaMOverSigmaM_vs_sigmaM->cd(2);
  if (data) hist2D_data->Draw("box");
  txt->DrawLatex(0.15,0.82,"Data");
  c_deltaMOverSigmaM_vs_sigmaM->cd(3);
  //hist2D_bkg->ProfileX()->Draw();
  hist2D_sig->Draw("box");
  txt->DrawLatex(0.15,0.82,"Signal");
  c_deltaMOverSigmaM_vs_sigmaM->cd(4);
  if (data) hist2D_data->ProfileX()->Draw();
  txt->DrawLatex(0.15,0.82,"Data (profile)");

  c_deltaMOverSigmaM_vs_sigmaM->SaveAs(outdir+"deltaMOverSigmaM_vs_sigmaM.gif");
  */

  /*
  float R_low = hist_bkg[2][13]->GetMean()/hist_bkg[1][13]->GetMean();
  float R_high = hist_bkg[2][13]->GetMean()/hist_bkg[3][13]->GetMean();

  cout << "<HT> (low) = " << hist_bkg[1][13]->GetMean() << endl;
  cout << "<HT> (sig) = " << hist_bkg[2][13]->GetMean() << endl;
  cout << "<HT> (high) = " << hist_bkg[3][13]->GetMean() << endl;
  cout << "R_low = " << R_low << endl;
  cout << "R_high = " << R_high << endl;
  */
}

void SetHistogramErrors(TH1* h) {

  float sf = h->Integral()/h->GetEntries();
  for (int ibin=0; ibin<h->GetNbinsX()+1; ibin++) {
    float err = sqrt(h->GetBinContent(ibin)/sf);
    h->SetBinError(ibin,err*sf);
  }

}
