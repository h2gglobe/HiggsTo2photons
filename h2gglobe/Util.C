#include "Util.h"
#include "LoopAll.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#define DEBUG 0

Util::Util() {
  loops = new LoopAll();
  loops->setUtilInstance(this);
#include "branchdef/newclonesarray.h"
}

void Util::SetTypeRun(int t, const char* name) {

  typerun = t;

  if (t == 1) {
    makeOutputTree = 1;
    outputFileName = TString(name);
    histFileName  = TString("");
  } else {
    makeOutputTree = 0;
    outputFileName = TString("");
    histFileName  = TString(name);
  }
  
  if (makeOutputTree) {
    cout << "CREATE " << outputFileName<<endl;
    outputFile=new TFile(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)
	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == 0 || typerun == 2) {
  	if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == 1) {
  	if(DEBUG) cout<<"typerun is " << typerun <<endl;

    hasoutputfile=1;
    outputTree = new TTree("event","Reduced tree"); 
    outputTreePar = new TTree("global_variables","Parameters");
    
    loops->PhotonAnalysisReducedOutputTree();
  } 
}

void Util::SetOutputNames(const char* name, const char* name2) {

  if (name2 != "")
    makeOutputTree = 1;
  histFileName  = TString(name);
  outputFileName = TString(name2);
}

void Util::AddFile(std::string name,int type) {
  if(DEBUG) cout << "Adding file:  " << name << " of type " << type << endl;

  files.push_back(name);
  datatype.push_back(type);

 
  //look in map for type as a key already
  map<int,int>::iterator it;
  if(DEBUG) cout << "map ok" << endl;
  it = type2HistVal.find(type);
  if(DEBUG) cout << "map checked" << endl;
  if (it == type2HistVal.end()) { 
    type2HistVal[type] = ntypes;
    if(DEBUG) cout << "container to add" << endl;
    ntypes++;
    if(DEBUG) cout << "container added" << endl;
  }
  nfiles++;
  if(DEBUG) cout << "successful addition" << endl;
}

void Util::LoopAndFillHistos(TString treename) {

  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal "<<endl;

  loops->InitReal(this,typerun);
  
  Files.resize(files.size());  
  Trees.resize(files.size());  
  TreesPar.resize(files.size());  

  std::vector<std::string>::iterator it;
  std::vector<TTree*>::iterator it_tree;
  std::vector<TTree*>::iterator it_treepar;
  std::vector<TFile*>::iterator it_file;

  it 	  	= files.begin();
  it_file 	= Files.begin();
  it_tree	= Trees.begin();
  it_treepar    = TreesPar.begin();  



  for (;it!=files.end()
       ;it_file++,it_tree++,it_treepar++,it++){
 
    cout<<"LoopAndFillHistos: opening " << i << " " << *it <<endl;
    this->current = i;
    *it_file = TFile::Open((*it).c_str());

    tot_events=1;
    sel_events=1;
    if(typerun == 1) { //this is a reduce job

      if(*it_file)
	*it_treepar=(TTree*) (*it_file)->Get("global_variables");

      if(*it_treepar) {
	TBranch        *b_tot_events;
	TBranch        *b_sel_events;
	(*it_treepar)->SetBranchAddress("tot_events",&tot_events, &b_tot_events);
	(*it_treepar)->SetBranchAddress("sel_events",&sel_events, &b_sel_events);
	b_tot_events->GetEntry(0);
	b_sel_events->GetEntry(0);
      } else {
	cout<<"REDUCE JOB, no global_variables tree ... the C-step job must have crashed ... SKIP THE FILE"<<endl;
	tot_events=0;
	sel_events=0;
      }
    }

    if(tot_events!=0) {

      if(*it_file)
	*it_tree=(TTree*) (*it_file)->Get(treename);

      loops->Init(typerun, *it_tree);
    }

    loops->Loop(this);

    if(tot_events != 0) {
      (*it_tree)->Delete("");
    }
   
    // EDIT - Cannot close the first file since it is in use after 
    // file 0 
    if(*it_file && i>0)
      (*it_file)->Close();
    
    i++;
  }
  //now close the first File
  if(Files[0]) Files[0]->Close();

  loops->TermReal(typerun);
}

void Util::WriteHist() {
    loops->myWritePlot(this);
}
