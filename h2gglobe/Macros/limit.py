#!/usr/bin/env python

from limit_func import *
import sys    # exit
import time   # time accounting
import getopt # command line parser

def main():
    method = "ada"
    print sys.argv
    try:
        # retrive command line options
        shortopts  = "M:T:D:P:A:"
        #longopts   = ["mass=","trees=","depth=","pruning=","ada="]
        opts, args = getopt.getopt( sys.argv[1:], shortopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        sys.exit(1)

    mass        = 0#DEFAULT_MASS
    tree        = 0#DEFAULT_TREES
    depth       = 0#DEFAULT_DEPTH
    prune       = "0"#DEFAULT_PRUNE
    ada       = 64#DEFAULT_NODES
    print opts
    print args 
    for o, a in opts:
        if o in ("-M", "--mass"):
            mass = int(a)
        elif o in ("-T", "--tree"):
            tree = int(a)
        elif o in ("-D", "--depth"):
            depth = int(a)
        elif o in ("-P", "--prune"):
            prune = str(a)
        elif o in ("-A", "--ada"):
            ada = float(a)
    run_limit(method,mass,tree,depth,prune,ada)


def run_limit(method,mass,tree,depth,prune,ada):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    setROOT()
    filename = "./TMVA_"+str(method)+"_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+str(prune)+"_"+str(ada)+".root"
    bdtvarname =  "BDT_"+str(method)+"_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+str(prune)+"_"+str(ada)
    trainName = "TrainTree"
    testName = "TestTree"
    ROOT.gROOT.cd()
    h_sig_train,h_bkg_train = fill_hist(filename,trainName,bdtvarname)

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

    r_list = [0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
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

    #f = open(filename+".txt", 'w')
    #f.write( "r\tCLs\n")
    #for r,CLs in zip(r_list,CLs_list):
    #    f.write(str(r)+"\t"+str(CLs)+"\n")
    #f.close()

    out_filename = "./Output_"+str(method)+"_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+str(prune)+"_"+str(ada)+".root"
    
    outfile = ROOT.TFile(out_filename,"RECREATE")
    out_name =  "BDT_"+str(method)+"_"+str(mass)+"_"+str(tree)+"_"+str(depth)+"_"+str(prune)+"_"+str(ada)
    ar_r   = array.array('d',r_list)
    ar_CLs = array.array('d',CLs_list)
    gr = ROOT.TGraph(len(r_list),ar_r,ar_CLs)
#
    f1 = ROOT.TF1("f1","pol4",r_list[0],r_list[-1]);
    gr.Fit("f1","R");

    gr.Write(out_name)
    f1.Write(out_name+"_fit")
    outfile.Close()

if __name__ == "__main__":
    main()
