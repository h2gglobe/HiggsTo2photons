import ROOT
import numpy,sys,math
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)

ROOT.gROOT.ProcessLine(".L quadInterpolate.C+")
from ROOT import quadInterpolate

from optparse import OptionParser

r=ROOT.TRandom3(0)
lumistring = "4.67 fb^{-1}"

def plainBin(hist):
	nb = hist.GetNbinsX()
	h2 = ROOT.TH1F(hist.GetName()+"new","",nb,0,nb)
	for i in range (1,nb+1):
		h2.SetBinContent(i,hist.GetBinContent(i))
		h2.SetBinError(i,hist.GetBinError(i))
	h2.GetXaxis().SetNdivisions(nb)
	return h2

def plotDistributions(mass,data,signals,bkg,errors):

	sigscale = 5.
	for i in range(1,len(signals)):
		signals[0].Add(signals[i])

	nbins = data.GetNbinsX()

	flatdata = plainBin(data)
	flatsignal  = plainBin(signals[0])
	flatbkg  = plainBin(bkg)

	fNew  = flatbkg.Clone()
	fNew2 = flatbkg.Clone()

	flatdata.SetMarkerStyle(20)
	flatdata.SetMarkerSize(1.0)
		
	fNew.SetLineColor(1)
	fNew.SetLineWidth(2)
	fNew.SetFillStyle(1001)
	fNew.SetFillColor(3)
	
	fNew2.SetFillColor(5)
	fNew2.SetLineColor(1)
	fNew2.SetLineWidth(2)
	fNew2.SetFillStyle(1001)
	fNewT = fNew.Clone()
	fNew2T = fNew2.Clone()
	fNewT.SetFillStyle(1001)
	fNew2T.SetFillStyle(1001)

	flatsignal.SetLineWidth(2)
	flatsignal.SetLineColor(ROOT.kRed)
	flatsignal.SetFillColor(ROOT.kRed-10)
	flatsignal.SetFillStyle(1001)
	flatsignal.Scale(sigscale)
		
	leg = ROOT.TLegend(0.56,0.56,0.88,0.88)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	leg.AddEntry(flatdata,"Data","PL")
	leg.AddEntry(flatsignal,"Higgs, m_{H}=%3.0f GeV (x5)"%(mass) ,"L")
	leg.AddEntry(flatbkg,"Bkg Model","L")
	leg.AddEntry(fNewT,"\pm 1\sigma","F")
	leg.AddEntry(fNew2T,"\pm 2\sigma","F")

	for b in range(1,nbins+1):
		additional = errors[b-1]
  		fNew.SetBinError(b,((fNew.GetBinError(b)**2)+((fNew.GetBinContent(b)*additional)**2))**0.5)
  		fNew2.SetBinError(b,2*(((fNew2.GetBinError(b)**2)+((fNew2.GetBinContent(b)*additional)**2))**0.5))
	
	c = ROOT.TCanvas("c","c",900,900)
	c.SetLogy()
	flatdata.GetXaxis().SetTitle("BDT bin number")
	flatdata.Draw("9")
	fNew2.Draw("9sameE2")
	fNew.Draw("9sameE2")
	flatbkg.Draw("9samehist")
	flatsignal.Draw("9samehist")
	flatdata.Draw("9sameP")

	leg.Draw()
	mytext = ROOT.TLatex()
	mytext.SetTextSize(0.03)
	mytext.SetNDC()
#	mytext.DrawLatex(0.28,0.8,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = %s}"%(lumistring))
	mytext.DrawLatex(0.1,0.92,"CMS preliminary,  #sqrt{s} = 7 TeV L = %s"%(lumistring))
	leg.Draw()
	c.SaveAs("m%3.1f.pdf"%mass)

def getBinningMass(mass):

	if mass >= 115.0 and mass <= 117.0: return "115"
	if mass >= 117.5 and mass <= 122.0: return "120"
	if mass >= 122.5 and mass <= 127.0: return "125"
	if mass >= 127.5 and mass <= 132.0: return "130"
	if mass >= 132.5 and mass <= 137.0: return "135"
	if mass >= 137.5 and mass <= 144.5: return "140"
	if mass >= 145.0 and mass <= 150.0: return "150"

def py_quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3):
	resL = quadInterpolate(-1*C,X1,X2,X3,Y1,Y2,Y3)
	resH = quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3)
	if math.isnan(resL) or math.isinf(resL) or  math.isnan(resH) or math.isinf(resL): return " - "
	if abs(resL - 1) < 0.00001 or abs(resL - 1) > 1: return " - "
	if abs(resH - 1) < 0.00001 or abs(resH - 1) > 1: return " - "
	return " %.8f/%.8f "%(resL,resH) 

def getBinContent(hist,b):
  
	res = hist.GetBinContent(b)
	if res==0: return 0.00000001
	else: return res

def getPoissonBinContent(hist,b,exp):
  
	res = exp*(hist.GetBinContent(b))
	if res==0: return 0.00000001
	else: return r.Poisson(res)

def writeCard(tfile,mass,scaleErr):

  print "Writing Datacard for mass -> ", mass
  outPut = open("mva-datacard_%3.1f.txt"%mass,"w")

  # Get All of the histograms we are going to use
  # Data ->
  dataHist = tfile.Get("th1f_data_grad_%3.1f_cat0"%mass)
  nBins    = dataHist.GetNbinsX()
  print "Number of Channels -> ", nBins
  # bkg model ->
  bkgHist  	= tfile.Get("th1f_bkg_grad_%3.1f_cat0"%mass)
  if options.Bias: bkgHistCorr   = tfile.Get("th1f_bkg_grad_%3.1f_cat0_fitsb_biascorr"%mass)
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

  if options.throwToy:
	print "Throwing toy dataset"
	for b in range(1,nBins+1): 
	  nd = r.Poisson(bkgHist.GetBinContent(b))
	  ns = 0
	  if options.expSig>0:
		print "Injecting %.f x SM"%expSig
		ns+=getPoissonBinContent(gghHist,b,options.expSig)
		ns+=getPoissonBinContent(vbfHist,b,options.expSig)
		ns+=getPoissonBinContent(wzhHist,b,options.expSig)
		ns+=getPoissonBinContent(tthHist,b,options.expSig)

	  outPut.write(" %.2f "%(nd+ns))

  else:
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

  backgroundContents = [bkgHist.GetBinContent(b) for b in range(1,nBins+1)]
  if options.Bias:
	print "Using Bkg Model Corrected for mass bias"
	backgroundContents = [bkgHistCorr.GetBinContent(b) for b in range(1,nBins+1)]

  for b in range(1,nBins+1): outPut.write(" %.12f   %.12f   %.12f   %.12f   %.12f "\
    %(getBinContent(gghHist,b),getBinContent(vbfHist,b),getBinContent(wzhHist,b),getBinContent(tthHist,b)\
     ,backgroundContents[b-1]))
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
  systematics = ["E_res","E_scale","idEff","r9Eff","kFactor","triggerEff","vtxEff"]
  
  if options.signalSys:
   print "Writing Systematics Part (coule be slow)"
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

  ## now for the David errors
  if options.Bias:
	print "Including Mass Bias nuisances (bin-to-bin stat error included)"

	# Input Signed Error Matrix from Fit 
	th2f_errmatrix = tfile.Get("fUncorrErr_%3.1f"%mass)
	for b in range(1,nBins+1):
           outPut.write("\nmassBias%d lnN"%b)
	   for q in range(1,nBins+1):
	   	f_errentry = th2f_errmatrix.GetBinContent(b,q)
		bkgC = backgroundContents[q-1]/sum(backgroundContents)
		bias_nuis  = py_quadInterpolate(1.,-1.,0.,1.,bkgC-f_errentry,bkgC,bkgC+f_errentry)
		outPut.write(" - - - - %s "%(bias_nuis))
        	
        outPut.write("\n")
	

  if options.B2B and (not options.Bias):
   # bkg bins will be gmN errors instead 
   for b in range(1,nBins+1):
        bkgScale = bkgHist.Integral()/bkgHist.GetEntries()
        outPut.write("\nbkg_stat%d gmN %d "%(b,int(backgroundContents[b-1]/bkgScale)))
	for q in range(1,nBins+1):
		if q==b: outPut.write(" - - - - %.8f "%bkgScale)
		else:    outPut.write(" - - - - - ")

  # Finally make a plot of what will go into the limit
  if options.makePlot: 
	print "Plotting Limit Setting inputs to m%3.1f.pdf"%mass
	if options.Bias : plotBKG = bkgHistCorr.Clone()
	else: plotBKG = bkgHist.Clone()
	plotDistributions(mass,dataHist.Clone() \
			,[gghHist.Clone(),vbfHist.Clone(),wzhHist.Clone(),tthHist.Clone()] \
			,plotBKG \
			,[scaleErr-1. for b in range(nBins)])

  outPut.write("\n")  # done
  outPut.close()
    
#################################################################################  

parser = OptionParser()
parser.add_option("-i","--input",dest="tfileName")
parser.add_option("","--noB2B",action="store_false",dest="B2B",default=True)
#parser.add_option("","--addBias",dest="biasFile",default=None)
parser.add_option("","--noBias",dest="Bias",default=True,action="store_false")
parser.add_option("","--noSignalSys",action="store_false",dest="signalSys",default=True)
parser.add_option("","--throwToy",action="store_true",dest="throwToy",default=False)
parser.add_option("","--expSig",dest="expSig",default=-1.,type="float")
parser.add_option("","--makePlot",dest="makePlot",default=False,action="store_true")
parser.add_option("-m","--mass",dest="singleMass",default=-1.,type="float")

(options,args)=parser.parse_args()
print "Creating Binned Datacards from workspace -> ", options.tfileName
if options.throwToy: print ("Throwing Toy dataset from BKG")
if options.expSig > 0: print ("(Also throwing signal SMx%f)"%options.expSig)

#if options.biasFile:
#	biasROOTFile = ROOT.TFile(options.biasFile)

genMasses     = [115,120,125,130,135,140,145,150]
#scalingErrors = [1.008,1.008,1.008,1.008,1.008,1.009,1.01,1.011] # Takes from P.Dauncey studies -> 7% window
#scalingErrors =  [1.007,1.007,1.006,1.008,1.007,1.008,1.009,1.01] # Takes from P.Dauncey studies -> 2% window
#scalingErrors = [1.013,1.013,1.012,1.012,1.014,1.015,1.016,1.016] # Takes from P.Dauncey studies -> 7% window (100-180)
#scalingErrors = [1.011,1.01,1.009,1.011,1.011,1.013,1.014,1.014] 	  # Takes from P.Dauncey studies -> 2% window (100-180)
scalingErrors = [1.00815,1.01024,1.01076,1.01197,1.0099,1.009,1.00928,1.01054 ] 	  # Takes from P.Dauncey studies -> 2% window (100-180) / MIT Preselection

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
print "Check the Errors Look Sensible -> plot saved to normErrors_%s"%options.tfileName
can.SaveAs("normErrors_%s.pdf"%options.tfileName)

# Now we can write the cards
#tfileName = sys.argv[1]
tfile = ROOT.TFile(options.tfileName)
if options.singleMass>0: evalMasses=[float(options.singleMass)]
for m in evalMasses: writeCard(tfile,m,normG.Eval(m))
print "Done Writing Cards"


