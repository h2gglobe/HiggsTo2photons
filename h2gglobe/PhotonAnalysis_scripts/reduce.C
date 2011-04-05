{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libRooFit.so");
  gSystem->Load("../libLoopAll.so");

  gBenchmark->Start("Reduction");
  Util* ut = new Util();

  ut->SetTypeRun(1, "GluGlu2H2GG140_reduced.root");
 // ut->AddFile("/vols/cms02/nw709/data.root",1);
  ut->AddFile("/vols/cms02/nw709/mc.root",1);
  //ut->AddFile("/vols/cms02/nw709/qcd_39.root",1);
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Reduction");
}

