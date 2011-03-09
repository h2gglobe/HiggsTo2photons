
{

  gSystem->Load("libRooFit.so");
  //using namespace RooFit;
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  gSystem->Load("../libLoopAll.so");

  cout << "Loaded libraries"<<endl;
  
  gBenchmark->Start("Analysis");
  Util* ut = new Util();

  cout << "Setting Analysis Type" << endl;
  ut->SetTypeRun(3, "stats.root");
  ut->AddFile("/vols/cms02/mjarvis/ntuples/h2gred/Run2010A.root",1);
  ut->AddFile("/vols/cms02/mjarvis/ntuples/h2gred/Run2010B.root",1);
  //ut->AddFileList("myFiles.txt");

  cout << "starting loop" <<endl;
  
  ut->LoopAndFillHistos();
  gBenchmark->Show("Anslysis");

//  ut->FitData();  
}

