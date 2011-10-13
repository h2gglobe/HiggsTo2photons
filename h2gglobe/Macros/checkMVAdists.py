import ROOT as r
r.gROOT.SetStyle("Plain")
r.gStyle.SetOptStat(0)
import sys
#r.gROOT.SetBatch(True)

print "Check MVA Signal Distributions -> Use as a quick check of the histograms going into the combination code (Run this over the interpolated workspace)"

rebin=False
if len(sys.argv) < 4: 
	sys.exit("checkMVAdists.py usage: python checkMVAdists.py FileName.root mass (eg 120) mode (eg grad) (option -rebin, rebins histogram to bin-number)")

fileName = sys.argv[1]
mass = float(sys.argv[2])
mode = str(sys.argv[3])

if len(sys.argv) ==5:
	if sys.argv[4]=="-rebin": rebin=True


f = r.TFile(fileName)
dat = f.Get("th1f_data_%s_%3.1f_cat0"%(mode,mass)).Clone()
bkg = f.Get("th1f_bkg_%s_%3.1f_cat0"%(mode,mass)).Clone()
sig = f.Get("th1f_sig_%s_%3.1f_cat0"%(mode,mass)).Clone()
datR = r.TH1F("datR","a",dat.GetNbinsX(),0,dat.GetNbinsX())
bkgR = r.TH1F("bkgR","b",bkg.GetNbinsX(),0,bkg.GetNbinsX())
sigR = r.TH1F("sigR","c",sig.GetNbinsX(),0,sig.GetNbinsX())

if rebin:
 for i in range(0,dat.GetNbinsX()+1):
  datR.SetBinContent(i,dat.GetBinContent(i));
  bkgR.SetBinContent(i,bkg.GetBinContent(i));
  sigR.SetBinContent(i,sig.GetBinContent(i));
  datR.SetBinError(i,dat.GetBinError(i));
  bkgR.SetBinError(i,bkg.GetBinError(i));
  sigR.SetBinError(i,sig.GetBinError(i));

 dat = datR.Clone()
 bkg = bkgR.Clone()
 sig = sigR.Clone()

sigInt = sig.Integral()
sig.Scale(bkg.Integral()/sig.Integral())

sig.SetLineColor(4)
sig.SetLineWidth(2)
bkg.SetFillColor(r.kOrange-3)
dat.SetMarkerStyle(23)
dat.SetMarkerSize(1.2)

sH = [sig.GetBinContent(i) for i in range(1,sig.GetNbinsX()+1)]
bH = [bkg.GetBinContent(i) for i in range(1,bkg.GetNbinsX()+1)]
dH = [dat.GetBinContent(i) for i in range(1,dat.GetNbinsX()+1)]
sH+=(bH)
sH+=(dH)

maxH = max(sH)
dat.SetMaximum(1.5*maxH)
dat.SetMinimum(0)

c = r.TCanvas("c","c",1000,1000)

leg = r.TLegend(0.15,0.6,0.65,0.8)
leg.SetFillColor(0)
leg.SetBorderSize(0)
dat.SetTitle("")
dat.GetYaxis().SetTitle("")
extra = ""
if mode == "grad": extra = "Gradient Boost"
if mode == "ada": extra = "Adaptive Boost"

if rebin: extra+=(" (bin number)")

dat.GetXaxis().SetTitle("BDT Output - %s"%extra)

leg.AddEntry(bkg,"Data Driven Background Model","F")
leg.AddEntry(sig,"Signal Model MC","L")
leg.AddEntry(dat,"Data","P")

bkgC = bkg.Clone()
bkgC.SetFillColor(2)
bkgC.SetFillStyle(4001)
bkgC.SetLineWidth(2)

dat.Draw()
bkg.Draw("samehist")
bkgC.Draw("samehistE2")
sig.Draw("samehist")
dat.Draw("sameP")
leg.Draw()


c.SaveAs("%s_mass_%3.1f.gif"%(mode,mass))
c.SaveAs("%s_mass_%3.1f.pdf"%(mode,mass))


d = r.TCanvas("d","d",1000,1000)
errh = bkg.Clone();
for i in range(1,errh.GetNbinsX()+1):
	if (dat.GetBinError(i) !=0):errh.SetBinContent(i,bkg.GetBinError(i)/dat.GetBinError(i))
	else: errh.SetBinContent(i,0)
	errh.SetBinError(i,0)

errh.GetYaxis().SetTitle("#sigma_{bkg}/#sigma_{obs}")
errh.SetLineWidth(2)
errh.SetMarkerStyle(21)
errh.SetMarkerSize(1.2)
errh.SetFillColor(0)
errh.SetMaximum(2)
errh.SetMinimum(0)
errh.Draw("hist")
errh.Draw("sameP")
d.SaveAs("%s_errs_mass_%d.gif"%(mode,mass))
