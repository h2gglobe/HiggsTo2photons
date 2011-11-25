import ROOT
import numpy,sys,math
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)

ROOT.gROOT.ProcessLine(".L quadInterpolate.C+")
from ROOT import quadInterpolate

def py_quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3):
	resL = quadInterpolate(-1*C,X1,X2,X3,Y1,Y2,Y3)
	resH = quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3)
	if math.isnan(resL) or math.isinf(resL) or  math.isnan(resH) or math.isinf(resL): return " - "
	if abs(resL - 1) < 0.001 or abs(resL - 1) > 1: return " - "
	if abs(resH - 1) < 0.001 or abs(resH - 1) > 1: return " - "
	#return " %.8f/%.8f "%(resL,resH) 
	#combination tool doesnt like kappas wrong way so lets symmetrize (take largest) but always give the UP direction
	diff = max(abs(1-resL),abs(1-resH))
	# now make direction, that of the up fluctuation
	sign = (resH-1)/abs(resH-1)
	return " %.8f "%(1+sign*diff)

def getBinContent(hist,b):
  
	res = hist.GetBinContent(b)
	if res==0: return 0.000001
	else: return res

def writeCard(tfile,mass,scaleErr):

  print "Writing Datacard for mass -> ", mass
  outPut = open("mva-datacard_%3.1f.txt"%mass,"w")

  # Get All of the histograms we are going to use
  # Data ->
  dataHist = tfile.Get("th1f_data_grad_%3.1f_cat0"%mass)
  nBins    = dataHist.GetNbinsX()
  # bkg model ->
  bkgHist  = tfile.Get("th1f_bkg_grad_%3.1f_cat0"%mass)
  # 4 signal channels ->
  gghHist  = tfile.Get("th1f_sig_grad_ggh_%3.1f_cat0"%mass)
  vbfHist  = tfile.Get("th1f_sig_grad_vbf_%3.1f_cat0"%mass)
  wzhHist  = tfile.Get("th1f_sig_grad_wzh_%3.1f_cat0"%mass)
  tthHist  = tfile.Get("th1f_sig_grad_tth_%3.1f_cat0"%mass)

 
  ###############################################################################
  # Write the basics
  outPut.write("Mva Binned Analysis DataCard (mH=%3.1f) \n"%mass)
  outPut.write("--------------------------------------------------------------\n")
  outPut.write("imax *\n")
  outPut.write("jmax *\n")
  outPut.write("kmax *\n")
  outPut.write("--------------------------------------------------------------\n")
  ## Now its the observation
  outPut.write("bin 	     ")
  for b in range(1,nBins+1): outPut.write(" cat%d "%b)
  outPut.write("\nobservation")
  for b in range(1,nBins+1): outPut.write(" %d "%dataHist.GetBinContent(b))
  outPut.write("\n--------------------------------------------------------------\n")
  ## Now we do the signal and background parts
  outPut.write("bin 	     ")
  for b in range(1,nBins+1): outPut.write(" cat%d  cat%d  cat%d  cat%d  cat%d "%(b,b,b,b,b))
  outPut.write("\nprocess    ")
  for b in range(1,nBins+1): outPut.write("  ggh    vbf    wzh    tth    bkg  ")
  outPut.write("\nprocess    ")
  for b in range(1,nBins+1): outPut.write("   0      0      0      0    1    ")
  outPut.write("\nrate       ")
  for b in range(1,nBins+1): outPut.write(" %.8f   %.8f   %.8f   %.8f   %.8f "%(getBinContent(gghHist,b),getBinContent(vbfHist,b),getBinContent(wzhHist,b),getBinContent(tthHist,b),bkgHist.GetBinContent(b)))
  outPut.write("\n--------------------------------------------------------------\n")

  # Some Globals #################################################################
  lumi 		= "1.045"
  QCDscale_ggH  = "0.918/1.125"
  PDF_gg_1      = "0.923/1.079"
  PDF_gg_2      = "0.915/1.085"
  QCDscale_qqH  = "0.997/1.005"
  PDF_qqbar_1   = "0.979/1.027"
  PDF_qqbar_2   = "0.958/1.042" 
  QCDscale_VH   = "0.982/1.018"
  QCDscale_ttH  = "0.905/1.036"
  ################################################################################

  # This next bit is for the signal systematics, first lets do the easy ones, lumi and theory
  outPut.write("\nlumi          lnN ")
  for b in range(1,nBins+1): outPut.write(" %s  %s  %s  %s  -  "%(lumi,lumi,lumi,lumi))
  outPut.write("\nQCDscale_ggH  lnN ")
  for b in range(1,nBins+1): outPut.write(" %s  -   -   -   -  "%(QCDscale_ggH))
  outPut.write("\nPDF_gg        lnN ")
  for b in range(1,nBins+1): outPut.write(" %s  -   -   %s  -  "%(PDF_gg_1,PDF_gg_2))
  outPut.write("\nQCDscale_qqH  lnN ")
  for b in range(1,nBins+1): outPut.write(" -   %s  -   -   -  "%(QCDscale_qqH))
  outPut.write("\nPDF_qqbar     lnN ")
  for b in range(1,nBins+1): outPut.write(" -   %s  %s  -   -  "%(PDF_qqbar_1,PDF_qqbar_2))
  outPut.write("\nQCDscale_VH   lnN ")
  for b in range(1,nBins+1): outPut.write(" -   -   %s  -   -  "%(QCDscale_VH))
  outPut.write("\nQCDscale_ttH  lnN ")
  for b in range(1,nBins+1): outPut.write(" -   -   -   %s  -  "%(QCDscale_ttH))

  outPut.write("\n")

  # Now is the very tedious part of the signal shape systematics, for each shape, simply do -/+ sigma
  print "Writing Systematics Part (coule be slow)"
  systematics = ["E_res","E_scale","idEff","r9Eff","kFactor","triggerEff","vtxEff"]
  
  for sys in systematics:

    gghHistU  = tfile.Get("th1f_sig_grad_ggh_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    vbfHistU  = tfile.Get("th1f_sig_grad_vbf_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    wzhHistU  = tfile.Get("th1f_sig_grad_wzh_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    tthHistU  = tfile.Get("th1f_sig_grad_tth_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    gghHistD  = tfile.Get("th1f_sig_grad_ggh_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    vbfHistD  = tfile.Get("th1f_sig_grad_vbf_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    wzhHistD  = tfile.Get("th1f_sig_grad_wzh_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    tthHistD  = tfile.Get("th1f_sig_grad_tth_%3.1f_cat0_%sDown01_sigma"%(mass,sys))

    outPut.write("\n%s lnN "%sys)
    for b in range(1,nBins+1): 
	 outPut.write(" %s %s %s %s - "%(\
				     py_quadInterpolate(1.,-3.,0.,3.,gghHistD.GetBinContent(b)  \
				        			  ,gghHist.GetBinContent(b)  \
                                        			  ,gghHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,vbfHistD.GetBinContent(b)  \
				        			  ,vbfHist.GetBinContent(b)  \
                                        			  ,vbfHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,wzhHistD.GetBinContent(b)  \
				        			  ,wzhHist.GetBinContent(b)  \
                                        			  ,wzhHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,tthHistD.GetBinContent(b)  \
				        			  ,tthHist.GetBinContent(b)  \
                                        			  ,tthHistU.GetBinContent(b)) \
 				    ))
  outPut.write("\n")
  # Finally the background errors, these are realtively simple
  outPut.write("\nbkg_norm lnN ")
  for b in range(1,nBins+1): outPut.write(" -   -   -   -  %.8f "%(scaleErr))

  # bkg bins will be gmN errors 
  bkgScale = bkgHist.Integral()/bkgHist.GetEntries()
  for b in range(1,nBins+1):
        outPut.write("\nbkg_stat%d gmN %d "%(b,int(bkgHist.GetBinContent(b)/bkgScale)))
	for q in range(1,nBins+1):
		if q==b: outPut.write(" - - - - %.8f "%bkgScale)
		else:    outPut.write(" - - - - - ")
  outPut.write("\n")  # done
  outPut.close()
  

                     
  
#################################################################################  
print "Creating Binned Datacards from workspace -> ", sys.argv[1]

genMasses     = [115,120,125,130,135,140,145,150]
scalingErrors = [1.00819,1.00764,1.00776,1.00794,1.00848,1.00917,1.0099,1.01143] # Takes from P.Dauncey studies
evalMasses    = numpy.arange(115,150.5,0.5)
normG = ROOT.TGraph(len(genMasses))

# Fill the errors graph
can = ROOT.TCanvas()
for i,ne in enumerate(scalingErrors):
  normG.SetPoint(i,genMasses[i],ne)
normG.SetMarkerStyle(20)
normG.GetXaxis().SetTitle("mH")
normG.GetYaxis().SetTitle("(N+dN)/N")
normG.Draw("ALP")
print "Check the Errors Look Sensible -> plot saved to normErrors_%s"%sys.argv[1]
can.SaveAs("normErrors_%s.pdf"%sys.argv[1])

# Now we can write the cards
tfileName = sys.argv[1]
tfile = ROOT.TFile(tfileName)
for m in evalMasses: writeCard(tfile,m,normG.Eval(m))
print "Done Writing Cards"


