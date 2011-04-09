{
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("libRooFit.so");
  gSystem->Load("libLoopAll.so");

  gBenchmark->Start("Analysis");
  Util* ut = new Util();

  //ut->SetTypeRun(2, "Run2010.root");
  //ut->AddFile("/vols/cms02/mjarvis/ntuples/h2gred/Run2010A.root", 1);
  //ut->AddFile("/vols/cms02/mjarvis/ntuples/h2gred/Run2010B.root", 1);
  
  ut->SetTypeRun(2, "QCD40.root");
  ut->AddFile("/vols/cms02/mjarvis/ntuples/h2gred/QCD40.root", 1);
   
  ut->LoopAndFillHistos();
  gBenchmark->Show("Analysis");

  ut->WriteHist();  
}

