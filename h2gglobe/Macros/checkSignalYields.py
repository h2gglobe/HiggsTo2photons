import ROOT
import sys,numpy

f = sys.argv[1]
if len(sys.argv)<2: 
	sys.exit("checkSignalYields.py usage:   python checkSignalYields.py NameOfFile.root (-interp for interpolated workspaces)")
doInterp=False
if len(sys.argv)>2:
	if sys.argv[2]=="-interp":
		doInterp=True

fil = ROOT.TFile(f)
print "Checking Signal / Background for - ", f
#m = [115,120,121,123,125,130,135,140,150]

if doInterp:
  m = numpy.arange(115,150.5,0.5)
  for mass in m :
	h = fil.Get("th1f_sig_ada_%3.1f_cat0"%mass)
	hB = fil.Get("th1f_bkg_ada_%3.1f_cat0"%mass)
	print "Mass-%3.1f Signal     (ada) -> "%mass ,h.Integral(), "Mass-%3.1f Background (ada) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()
  print 
  for mass in m :
	h = fil.Get("th1f_sig_grad_%3.1f_cat0"%mass)
	hB = fil.Get("th1f_bkg_grad_%3.1f_cat0"%mass)
	print "Mass-%3.1f Signal     (grad) -> "%mass ,h.Integral(), "Mass-%3.1f Background (grad) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()

else:
  m = [115,120,125,130,135,140,150]
  for mass in m :
	h = fil.Get("th1f_sig_ada_%d.0_%d.0_cat0"%(mass,mass))
	hB = fil.Get("th1f_bkg_ada_%d.0_cat0"%(mass))
	print "Mass-%d Signal     (ada) -> "%mass ,h.Integral(), "Mass-%d Background (ada) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()
  print
  for mass in m :
	h = fil.Get("th1f_sig_grad_%d.0_%d.0_cat0"%(mass,mass))
	hB = fil.Get("th1f_bkg_grad_%d.0_cat0"%(mass))
	print "Mass-%d Signal     (grad) -> "%mass ,h.Integral(), "Mass-%d Background (grad) -> "%mass ,hB.Integral(), " S/B ",h.Integral()/hB.Integral()
