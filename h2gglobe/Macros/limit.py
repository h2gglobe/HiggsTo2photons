#!/usr/bin/env python

from limit_func import *
import sys    # exit
import time   # time accounting
import getopt # command line parser

def main():
    try:
        # retrive command line options
        shortopts  = "M:T:D:P:N"
        longopts   = ["mass=","trees=","depth=","pruning=","nodes="]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        sys.exit(1)

    mass        = 0#DEFAULT_MASS
    tree        = 0#DEFAULT_TREES
    depth       = 0#DEFAULT_DEPTH
    prune       = "0"#DEFAULT_PRUNE
    node       = 64#DEFAULT_NODES
    print opts
    for o, a in opts:
        if o in ("-M", "--mass"):
            mass = int(a)
        elif o in ("-T", "--tree"):
            tree = int(a)
        elif o in ("-D", "--depth"):
            depth = int(a)
        elif o in ("-P", "--prune"):
            prune = str(a)
        elif o in ("-N", "--nodes"):
            node = int(a)
    run_limit(mass,tree,depth,prune,node)


def run_limit(mass,tree,depth,prune,node):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    #mass_list = [120] 
    #tree_list = [128,256,512]
    #depth_list= [4,8,16,32]
    #prune_list = ["NoPruning"]#,"CostComplexity","ExpectedError"]
    #node_list = [4,8,16,32,64,128,1024]
    setROOT()
    #for mass in mass_list : 
    #    for tree in tree_list : 
    #        for depth in depth_list :
    #            for prune in prune_list : 
    #                for node in node_list : 
    #                    call(["qsub", "-q", "hepshort.q", "subTMVAToBatch.sh", str(pwd) , str(mass), str(tree), str(depth), str(prune), str(node)])
    filename = "./TMVA_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+prune+"_"+str(node)+".root"
    bdtvarname = "BDT_grad_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+prune#TODO correct this to the actual variable name
    #bdtvarname = "./TMVA_120_128_16_NoPruning_1024.root"
    trainName = "TrainTree"
    testName = "TestTree"
    ROOT.gROOT.cd()
    h_sig_train,h_bkg_train = fill_hist(filename,trainName,bdtvarname)
    print "\ts :", h_sig_train.Integral()
    print "\tb :", h_bkg_train.Integral()

    arBins = opt_binning(h_bkg_train,50,False,True)
    #def opt_binning(hb,nTargetBins,revise_target,use_n_target):

    # rebin histograms and output plots to esure everything is correct
    h_bkg_train_opt = h_bkg_train.Rebin(len(arBins)-1,"h_bkg_train_opt",arBins)
    c = ROOT.TCanvas("c","c",800,600)
    c.cd()
    h_bkg_train.Draw()
    c.SaveAs("./raw.pdf")
    h_bkg_train_opt.Draw()
    c.SaveAs("./opt.pdf")

    h_sig_test,h_bkg_test = fill_hist(filename,testName,bdtvarname) 

    h_bkg_test_opt = h_bkg_test.Rebin(len(arBins)-1,"h_bkg_test_opt",arBins)
    h_sig_test_opt = h_sig_test.Rebin(len(arBins)-1,"h_sig_test_opt",arBins)

    h_data = h_bkg_test_opt.Clone("h_data")

    c.cd()
    h_bkg_test_opt.Draw()
    c.SaveAs("./test_bkg.pdf")
    h_sig_test_opt.Draw()
    c.SaveAs("./test_sig.pdf")

    h_sig_test,h_bkg_test = fill_hist(filename,testName,bdtvarname) 

    r_list = [0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
    CLs_list = []

    lim = ROOT.TLimit()
    for r in r_list :
        h_sig_test_opt.Scale(r)  
        print "Calculating CL_s at r =", r
        print "\ts :", h_sig_test_opt.Integral()
        print "\tb :", h_bkg_test_opt.Integral()
        cl = lim.ComputeLimit(h_sig_test_opt,h_bkg_test_opt,h_data) 
        CLs_list.append(1.-cl.GetExpectedCLs_b())  
        h_sig_test_opt.Scale(1/r)  

    f = open(filename+".txt", 'w')
    f.write( "r\tCLs\n")
    for r,CLs in zip(r_list,CLs_list):
        f.write(str(r)+"\t"+str(CLs)+"\n")
    f.close()




if __name__ == "__main__":
    main()
