/* Background Shape Systematic for use with h2gglobe MvaAnalysis 
Original Author - D. Futyan
*/

#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iomanip>
#include <string>


void bkgShapeSyst_mass(TFile *f_in,std::string mode, double mass, int cat,bool equalBinWidths=false) {


  TH1* hist_data[4];
  TH1* hist_data_norm[4];

  TString title, var;
  if (equalBinWidths) {
    title=Form("BDT output bin number (%s)",mode.c_str());
    var=Form("bkgSyst_binNum_%s_syst",mode.c_str());
  } else {
    title=Form("BDT output (%s)",mode.c_str());
    var=Form("bkgSyst_%s",mode.c_str());
  }

//  TFile *f_out = TFile::Open(Form(var+"_%d.root","recreate");

  f_in->cd();

  if (!equalBinWidths) {
    hist_data[1] = (TH1*)f_in->Get(Form("th1f_bkg_1low_%s_%3.1f_cat%d",mode.c_str(),mass,cat))->Clone();
    hist_data[2] = (TH1*)f_in->Get(Form("th1f_bkg_%s_%3.1f_cat%d",mode.c_str(),mass,cat))->Clone();
    hist_data[3] = (TH1*)f_in->Get(Form("th1f_bkg_1high_%s_%3.1f_cat%d",mode.c_str(),mass,cat))->Clone();
  //  hist_data[2] = (TH1*)th1f_bkg_ada_120_cat0->Clone();
  //  hist_data[3] = (TH1*)th1f_bkg_high_ada_120_cat0->Clone();
  } else {
    TH1 *th1f_bkg_ada_120_cat0 = (TH1*) f_in->Get(Form("th1f_bkg_%s_%3.1f_cat%d",mode.c_str(),mass,cat));
    TH1 *th1f_bkg_high_ada_120_cat0 = (TH1*) f_in->Get(Form("th1f_bkg_1high_%s_%3.1f_cat%d",mode.c_str(),mass,cat));
    TH1 *th1f_bkg_low_ada_120_cat0 = (TH1*) f_in->Get(Form("th1f_bkg_1low_%s_%3.1f_cat%d",mode.c_str(),mass,cat));
    int nbins=th1f_bkg_ada_120_cat0->GetNbinsX();
    hist_data[1] = new TH1F(Form("hist_data_low_%s_%3.1f_cat%d",mode.c_str(),mass,cat),"hist_data_low",nbins,0,float(nbins));
    hist_data[2] = new TH1F(Form("hist_data_sig_%s_%3.1f_cat%d",mode.c_str(),mass,cat),"hist_data_sig",nbins,0,float(nbins));
    hist_data[3] = new TH1F(Form("hist_data_high_%s_%3.1f_cat%d",mode.c_str(),mass,cat),"hist_data_high",nbins,0,float(nbins));
    for (int ibin=0; ibin<nbins+1; ibin++) {
      hist_data[1]->SetBinContent(ibin,th1f_bkg_low_ada_120_cat0->GetBinContent(ibin));
      hist_data[1]->SetBinError(ibin,th1f_bkg_low_ada_120_cat0->GetBinError(ibin));
      hist_data[2]->SetBinContent(ibin,th1f_bkg_ada_120_cat0->GetBinContent(ibin));
      hist_data[2]->SetBinError(ibin,th1f_bkg_ada_120_cat0->GetBinError(ibin));
      hist_data[3]->SetBinContent(ibin,th1f_bkg_high_ada_120_cat0->GetBinContent(ibin));
      hist_data[3]->SetBinError(ibin,th1f_bkg_high_ada_120_cat0->GetBinError(ibin));
    }
  }

  for (int i=1; i<4; i++) {
    hist_data[i]->SetMarkerStyle(20);
    hist_data[i]->SetMarkerSize(.8);
    hist_data[i]->GetXaxis()->SetTitle(title);
    hist_data[i]->GetXaxis()->SetTitleSize(0.04);
    hist_data[i]->GetYaxis()->SetTitle("");
  }

  for (int i=1; i<4; i++) {
    float sf=hist_data[i]->Integral();
    hist_data_norm[i] = (TH1*)hist_data[i]->Clone();
    hist_data_norm[i]->Scale(1./sf);
  }

  TH1* hist_data_highMinusLow = (TH1*)hist_data_norm[1]->Clone();
  hist_data_highMinusLow->Add(hist_data_norm[3],-1.);
  hist_data_highMinusLow->Scale(0.5);

  TH1* th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma = (TH1*)hist_data_norm[2]->Clone();
  TH1 *th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma = (TH1*)hist_data_norm[2]->Clone();
  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->SetName(Form("th1f_bkg_%s_%3.1f_cat%d_bkgShapeUp01_sigma",mode.c_str(),mass,cat));
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->SetName(Form("th1f_bkg_%s_%3.1f_cat%d_bkgShapeDown01_sigma",mode.c_str(),mass,cat));

  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->Add(hist_data_highMinusLow);
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->Add(hist_data_highMinusLow,-1);

  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->Scale(hist_data[2]->Integral());
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->Scale(hist_data[2]->Integral());

  gROOT->SetBatch(true);

  TCanvas *canvas = new TCanvas(Form("c_%s%3.1f",mode.c_str(),mass)+var,var);
  canvas->SetFillColor(0);

  hist_data[2]->Draw();
  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->SetMarkerColor(2);
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->SetMarkerColor(kGreen-2);
  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->SetLineColor(2);
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->SetLineColor(kGreen-2);
  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->Draw("same");
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->Draw("same");

  TLegend *leg = new TLegend(.15,.75,.47,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_data[2],"Background model: nominal","L");
  leg->AddEntry(th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma,"Background model: up 1 sigma","L");
  leg->AddEntry(th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma,"Background model: down 1 sigma","L");
  leg->Draw();

  canvas->SaveAs(Form("bkgSystematic_%s_%3.1f_cat%d.pdf",mode.c_str(),mass,cat));
  f_in->cd();
  th1f_bkg_ada_120_cat0_bkgShapeUp01_sigma->Write();
  th1f_bkg_ada_120_cat0_bkgShapeDown01_sigma->Write();

}

void bkgShapeSyst(std::string FileName, int ncat=1,bool equalBinWidths=false){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  TFile *f_in =(TFile*) TFile::Open(FileName.c_str(),"UPDATE");

//  const int nmasses = 9;
//  int masses[nmasses] = {115,120,121,123,125,130,135,140,150};
//  int masses[nmasses] = 
  //for (int m=0;m<nmasses;m++){
  for (double m=115.;m<=150.;m+=0.5){
    for (int cat=0;cat<ncat;cat++){
       bkgShapeSyst_mass(f_in,"ada",m,cat,equalBinWidths);
       bkgShapeSyst_mass(f_in,"grad",m,cat,equalBinWidths);
    }
  }


}
