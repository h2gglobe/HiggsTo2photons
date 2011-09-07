import ROOT
import sys
import numpy
from theory_sm import *

ROOT.gROOT.SetStyle("Plain")
f = ROOT.TFile(sys.argv[1])
lumi = 1658.
dataName ="th1f_data_mass"

print "File - ", sys.argv[1]
printLine = "Data:      "
Sum = 0
for i in range(8):
  h = f.Get(dataName+("_cat%d"%i))
  Sum+=h.Integral()
  printLine+="%3.0f"%(h.Integral())+" "
printLine+="tot=%d"%Sum

print printLine

efficiency=ROOT.TGraphAsymmErrors()
central=ROOT.TGraphAsymmErrors()
efficiencyup=ROOT.TGraphAsymmErrors()
efficiencydn=ROOT.TGraphAsymmErrors()
centralsmooth=ROOT.TGraphAsymmErrors()
MG=ROOT.TMultiGraph()
can = ROOT.TCanvas()
systematics = ["vtxEff","idEff","E_scale","E_res","triggerEff"]
# central Values
Masses = numpy.arange(110,151,1.)
#Masses=[110,120,130,140]
fitstring = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x"
cenfunc = ROOT.TF1("cenfunc",fitstring,109.75,140.25)
upfunc = ROOT.TF1("upfunc",fitstring,109.75,140.25)
dnfunc = ROOT.TF1("dnfunc",fitstring,109.75,140.25)


for point,M in enumerate(Masses):
  printLine = "Sigl M%d: "%M
  Sum = 0
  for i in range(8):
    if int(M)==M:
     h = f.Get("th1f_sig_mass_m%d_cat%d"%(int(M),i))
    else:
     h = f.Get("th1f_sig_mass_m%.1f_cat%d"%(M,i))
    Sum+=h.Integral()
    printLine+="%3.5f"%h.Integral()+" "
  printLine+="tot=%3.5f"%Sum
  
  sm =1.0
  for j,mm in enumerate(allMasses): 
	if mm==M: 
		sm = xSec[j]*br[j]
  effAcc = Sum/(sm*lumi)
  centralsmooth.SetPoint(point,M,effAcc)
  central.SetPoint(point,M,effAcc)
  efficiency.SetPoint(point,M,effAcc)
  delUp = 0
  delDown = 0
  for s in systematics:
   syssumup=0
   syssumdn=0
   for i in range(8):
    if int(M)==M:
     print M, s
     hup = f.Get("th1f_sig_mass_m%d_cat%d_%sUp01_sigma"%(int(M),i,s))
     hdn = f.Get("th1f_sig_mass_m%d_cat%d_%sDown01_sigma"%(int(M),i,s))
    else:
     print "M,s,th1f_sig_mass_m%.1f_cat%d_%sUp01_sigma"%(M,i,s)
     hup = f.Get("th1f_sig_mass_m%.1f_cat%d_%sUp01_sigma"%(M,i,s))
     hdn = f.Get("th1f_sig_mass_m%.1f_cat%d_%sDown01_sigma"%(M,i,s))
    syssumup+=hup.Integral()
    syssumdn+=hdn.Integral()
   delUp+=abs(syssumup-Sum)/3
   delDown+=abs(syssumdn-Sum)/3
   
  delUp=(delUp**0.5)/(sm*lumi)
  delDown=(delDown**0.5)/(sm*lumi)
  efficiencyup.SetPoint(point,M,delUp)
  efficiencydn.SetPoint(point,M,delDown)
  centralsmooth.SetPointError(point,0,0,0,0)
  efficiency.SetPointError(point,0,0,delDown,delUp)
#  efficiency.SetPointError(point,0,0,0,0)

  print printLine

centralsmooth.Fit(cenfunc,"R,0,EX0","")
efficiencyup.Fit(upfunc,"R,0,EX0","")
efficiencydn.Fit(dnfunc,"R,0,EX0","")

for point,M in enumerate(Masses):
 central.SetPoint(point,M,cenfunc.Eval(M))
 efficiency.SetPoint(point,M,cenfunc.Eval(M))
# efficiency.SetPointError(point,0,0,dnfunc.Eval(M),upfunc.Eval(M))

efficiency.SetFillColor(9)
central.SetLineWidth(2)
MG.Add(efficiency)
MG.Add(central)
MG.Draw("AL3")
MG.GetXaxis().SetTitle("m_{H}(GeV/c^{2})")
MG.GetYaxis().SetTitle("Efficiency #times Acceptance - %")
MG.GetYaxis().SetRangeUser(33,49)

leg=ROOT.TLegend(0.46,0.16,0.79,0.39)
#leg=ROOT.TLegend(0.46,0.62,0.79,0.86)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(central,"Higgs Signal #varepsilon #times Acc","L")
leg.AddEntry(efficiency,"#pm 1 #sigma syst. error","F")

mytext = ROOT.TLatex()
mytext.SetTextSize(0.04)
mytext.SetNDC()
mytext.DrawLatex(0.15,0.8,"CMS Simulation")
leg.Draw()
can.SaveAs("effAcc_vs_mass.pdf")

  

	
