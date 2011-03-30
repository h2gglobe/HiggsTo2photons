/**************************************************
// Utility to read cfg file dump into globe ntuples
//
// Matteo Sani, UCSD
// 30/03/2011
//
**************************************************/

read_global_variables(const char* filename) {

  std::string* jobmaker = new std::string();
  std::vector<std::string>* parameters = new std::vector<std::string>();
  Int_t version, tot_events, sel_events;
  
  TFile* file = TFile::Open(filename);
  TTree* tree = (TTree*)file->Get("global_variables");

  tree->SetBranchAddress("jobmaker", &jobmaker);
  tree->SetBranchAddress("parameters", &parameters);
  tree->SetBranchAddress("version", &version);
  tree->SetBranchAddress("tot_events", &tot_events);
  tree->SetBranchAddress("sel_events", &sel_events);

  tree->GetEntry(0);

  std::cout << "H2GAnalyzer Version: " << version << std::endl;  
  std::cout << "Processed events: " << tot_events << std::endl;  
  std::cout << "Selected events : " << sel_events << std::endl;  
  std::cout << *jobmaker << std::endl;
  std::cout << "CONFIGURATION: " << std::endl;
  for (unsigned int i=0; i<parameters->size(); i++) {
    std::cout << (*parameters)[i] << std::endl;
  }
}
