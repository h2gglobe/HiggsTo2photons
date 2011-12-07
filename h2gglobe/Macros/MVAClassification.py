
#!/usr/bin/env python
# @(#)root/tmva $Id: MVAClassification.py,v 1.1.4.4 2011/12/05 09:36:41 mjarvis Exp $
# ------------------------------------------------------------------------------
# based on TMVA Python script: TMVAClassification.py
# ------------------------------------------------------------------------------

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser

# --------------------------------------------

# Default settings for command line arguments
DEFAULT_OUTFNAME = "TMVA"
DEFAULT_INFNAME  = "TMVA_input_CMS-HGG_4686pb.root"
DEFAULT_TREESIG  = "sig"
DEFAULT_TREEBKG  = "bkg"
DEFAULT_METHODS  = "BDT"
DEFAULT_BACKGROUND = 1 # 1 uses the sidebands, 2 uses all
DEFAULT_MASS     = 123
DEFAULT_CAT      = -1
DEFAULT_WIDTH    = 0.02

# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: '%s')" % DEFAULT_METHODS  
    print "  -M | --mass       : gives Higgs Mass to train (default: '%i')" % DEFAULT_MASS  
    print "  -C | --cat        : Diphoton selection cat to run over (default: '%i')" % DEFAULT_MASS  
    print "  -B | --background : define which background samples to use (default: '%i' [sidebands])" % DEFAULT_BACKGROUND  
    print "  -i | --inputfile  : name of input ROOT file (default: '%s')" % DEFAULT_INFNAME
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" % (DEFAULT_TREESIG, DEFAULT_TREEBKG)
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "m:M:C:B:i:t:o:vh?"
        opts, args = getopt.getopt( sys.argv[1:], shortopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    infname     = DEFAULT_INFNAME
    methods     = DEFAULT_METHODS
    mass        = DEFAULT_MASS
    cat         = DEFAULT_CAT
    outfname    = DEFAULT_OUTFNAME
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    bkg_method  = DEFAULT_BACKGROUND
    width       = DEFAULT_WIDTH
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-M", "--mass"):
            mass = int(a)
        elif o in ("-C", "--cat"):
            cat = int(a)
        elif o in ("-B", "--background"):
            bkg_method = int(a)
        elif o in ("-i", "--inputfile"):
            infname = a
        elif o in ("-o", "--outputfile"):
            outfname = a
        elif o in ("-t", "--inputtrees"):
            a.strip()
            trees = a.rsplit( ' ' )
            trees.sort()
            trees.reverse()
            if len(trees)-trees.count('') != 2:
                print "ERROR: need to give two trees (each one for signal and background)"
                print trees
                sys.exit(1)
            treeNameSig = trees[0]
            treeNameBkg = trees[1]
        elif o in ("-v", "--verbose"):
            verbose = True

    if (width == 0.02) : width_str = "_2pt"
    elif (width == 0.07) : width_str = "_7pt"
    mass_str    = "_"+str(mass)
    cat_str    = "_"+str(cat)
    if cat<0:
        cat_str    = "_all"
    outfname    = outfname+mass_str+cat_str+".root"
    #treeNameSig = treeNameSig + mass_str 
    #treeNameBkg = treeNameBkg + mass_str  

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)." 
        sys.exit(1)
    
    # Logon not automatically loaded through PyROOT (logon loads TMVA library)
    # load also GUI
    gROOT.SetMacroPath( "./" )
    #gROOT.Macro       ( "./TMVAlogon.C" )    
    #gROOT.LoadMacro   ( "./TMVAGui.C" )
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( outfname, 'RECREATE' )

    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for
    # more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = TMVA.Factory( "TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification")

    # Set verbosity
    factory.SetVerbose( verbose )
    
    # Define the input variables that shall be used for the classifier training
    factory.AddVariable( "pho1_ptOverM", "P_{T}^{lead} / M_{H}", "", 'F' );
    factory.AddVariable( "pho2_ptOverM", "P_{T}^{sublead} / M_{H}", "", 'F' );
    factory.AddVariable( "pho1_eta","#eta^{lead}", "", 'F' );
    factory.AddVariable( "pho2_eta","#eta^{sublead}", "", 'F' );
    factory.AddVariable( "d_phi","#Delta #phi", "rad", 'F' );#should this be cos delta Phi
    factory.AddVariable( "H_ptOverM","P_{T}^{Higgs}/M_{H}", "", 'F' );
    factory.AddVariable( "H_eta","#eta^{Higgs}", "", 'F' );

    factory.AddVariable( "sigmaMOverM","#sigmaM_{cor} / M",  'F' )
    factory.AddVariable( "sigmaMOverM_wrongVtx","#sigmaM_{wrong} / M",  'F' )

    factory.AddVariable( "vtx_prob","P_{vertex}", "", 'F' );

    #factory.AddVariable( "cos_theta_star","cos(#theta)*", "", 'F' );
    #factory.AddVariable( "max_eta","max(#eta^{lead},#eta^{sub.})", "", 'F' );
    #factory.AddVariable( "min_r9","min(r9^{lead},r9^{sub.})", "", 'F' );

    factory.AddVariable( "deltaMOverM","#DeltaM / M_{Hypth}.",  'F' )
        
    input = TFile.Open( infname )

    # Get the signal and background trees for training
    signal_train      = input.Get( treeNameSig+"_train"+mass_str+".0")
    signal_test      = input.Get( treeNameSig+"_test"+mass_str+".0")

    background_train  = input.Get( treeNameBkg+"_train"+width_str+mass_str+".0")
    background_test  = input.Get( treeNameBkg+"_test"+width_str+mass_str+".0")

    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    # ====== register trees ====================================================
    factory.AddSignalTree    ( signal_train,    signalWeight     ,"train")
    factory.AddBackgroundTree( background_train,backgroundWeight ,"train")
    factory.AddSignalTree    ( signal_test,     signalWeight     ,"test")
    factory.AddBackgroundTree( background_test, backgroundWeight ,"test")
            
    # Set individual event weights (the variables must exist in the original
    # TTree)
    factory.SetBackgroundWeightExpression( "wt" )
    factory.SetSignalWeightExpression( "wt" )

    # Apply additional cuts on the signal and background sample. 
    # example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )
    mycut = TCut( "fabs(deltaMOverM)<="+str(width))#
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    factory.PrepareTrainingAndTestTree( mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:NormMode=NumEvents:!V")
    # Boosted Decision Trees
    # NEW PARAMETERS
    factory.BookMethod( TMVA.Types.kBDT, "BDT_ada" +mass_str+cat_str,"!H:!V:NTrees=400:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.05:SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning")
    factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+mass_str+cat_str,"!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=50:NNodesMax=5") 
    #test

    # --------------------------------------------------------------------------------------------------
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

    # Train MVAs
    #factory.OptimizeAllMethods()
    factory.TrainAllMethods()
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    #gROOT.ProcessLine( "TMVAGui(\"%s\")" % outfname )
    
    # keep the ROOT thread running
    #gApplication.Run() 

# ----------------------------------------------------------

if __name__ == "__main__":
    main()

