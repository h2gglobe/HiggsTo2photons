import ROOT
#from readJSON import ReadJSON
import sys
name = sys.argv[1]

ROOT.gSystem.Load("libRooFit.so");
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

#jsonFile = "fileList.JSON"
#json = ReadJSON(jsonFile)

ROOT.gBenchmark.Start("Analysis");
ut = ROOT.Util();

ut.SetTypeRun(2, 'hist/'+(name.split('/'))[-1]+'_hist.root');
ut.AddFile(name+'.root',1);
  
ut.LoopAndFillHistos();
ROOT.gBenchmark.Show("Analysis");

ut.WriteHist();  


