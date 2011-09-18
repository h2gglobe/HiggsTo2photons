import ROOT
import sys

grepmode=False
if len(sys.argv) < 2: 
	sys.exit("grabAllRooPlots.py usage:  python grabAllRooPlots.py FileName.root (optional -grep=expr , only plots containing expr in name are used)")

fileName = sys.argv[1]

if len(sys.argv) == 3:
  if "-grep=" in sys.argv[2]:
	grepmode = True
	grep	 = ((sys.argv[2]).split("="))[-1]

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetStyle("Plain")
F = ROOT.TFile(fileName)
keys = F.GetListOfKeys()
plots = []

for K in keys:
	obj = K.ReadObj()
	if grepmode:
	  print "grepping ", grep
	  if "plot" in obj.GetName() and grep in obj.GetName(): plots.append(obj.Clone())
	else:
	  if "plot" in obj.GetName(): plots.append(obj.Clone())

can = ROOT.TCanvas("c","Fit Plots",2000,2000)
nPlots = len(plots)
sqrtnPlots = int(nPlots**0.5)+1
can.Divide(sqrtnPlots,sqrtnPlots)

for i,P in enumerate(plots):
  can.cd(i+1)
  P.DrawClonePad()

can.SaveAs("allThePlots.eps")
