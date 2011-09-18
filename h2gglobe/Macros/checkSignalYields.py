import ROOT
import sys

f = sys.argv[1]
if len(sys.argv)!=2: 
	sys.exit("checkSignalYields.py usage:   python checkSignalYields.py NameOfFile.root")

fil = ROOT.TFile(f)
print "Checking Signal / Background for - ", f
m = [115,120,121,123,125,130,135,140,150]
for mass in m :
	h = fil.Get("th1f_sig_ada_%d_cat0"%mass)
	hB = fil.Get("th1f_bkg_ada_%d_cat0"%mass)
	print "Mass-%d Signal     (ada) -> "%mass ,h.Integral(), "Mass-%d Background (ada) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()
print
for mass in m :
	h = fil.Get("th1f_sig_grad_%d_cat0"%mass)
	hB = fil.Get("th1f_bkg_grad_%d_cat0"%mass)
	print "Mass-%d Signal     (grad) -> "%mass ,h.Integral(), "Mass-%d Background (grad) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()
