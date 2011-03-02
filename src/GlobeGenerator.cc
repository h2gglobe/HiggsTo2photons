#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenerator.h"

#include <iostream>

GlobeGenerator::GlobeGenerator(const edm::ParameterSet& iConfig) {
  
  generatorColl = iConfig.getParameter<edm::InputTag>("GeneratorColl");
  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  edm::ParameterSet psetGenerator = iConfig.getParameter<edm::ParameterSet>("GeneratorCuts");
  
  // get cut thresholds
  etCut_ =  psetGenerator.getParameter<double>("EtCut");
}

void GlobeGenerator::defineBranch(TTree* tree) {

  // think about changing branch names for duplicate collections
  gen_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
  
  tree->Branch("gen_n", &gen_n, "gen_n/I");
  
  tree->Branch("gen_p4", "TClonesArray", &gen_p4, 32000, 0);
  
  tree->Branch("gen_status", gen_status, "gen_status[gen_n]/I");
  tree->Branch("gen_pdgid", gen_pdgid, "gen_pdgid[gen_n]/I");
  tree->Branch("gen_mother", gen_mother, "gen_mother[gen_n]/I");
}

bool GlobeGenerator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  // take collections
  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByLabel(generatorColl, evt);
  myGenEvent = evt->GetEvent();

  gen_p4->Clear();
  
  gen_n = 0;
  
  if(debug_level > 99)
    std::cout << "#"
	      << "\t pdgid"
	      << "\t status"
	      << "\t mother"
	      << "\t px \t py \t pz \t E"
	      << "\t vx \t vy \t vz \t t"
	      << std::endl;

  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 
        
    if (gen_n >= MAX_GENERATOR) {
      std::cout << "GlobeGenerator: WARNING TOO MANY Generator PARTICLES (allowed " << MAX_GENERATOR << ")" << std::endl;
      break;
    }

    new ((*gen_p4)[gen_n]) TLorentzVector();
    ((TLorentzVector *)gen_p4->At(gen_n))->SetXYZT((*it)->momentum().px(), (*it)->momentum().py(),
                                                   (*it)->momentum().pz(), (*it)->momentum().e());
    gen_pdgid[gen_n] = (*it)->pdg_id();
    gen_status[gen_n] = (*it)->status();
    gen_mother[gen_n] = mother(myGenEvent, *it);
    
    if(debug_level > 99) {
      std::cout << gen_n
		<< "\t" << gen_pdgid[ gen_n]
		<< "\t" << gen_status[gen_n]
		<< "\t" << gen_mother[gen_n]
		<< "\t" << (*it)->momentum().px()
		<< "\t" << (*it)->momentum().py()
		<< "\t" << (*it)->momentum().pz()
		<< "\t" << (*it)->momentum().e();
      if((*it)->production_vertex())
	std::cout << "\t" << (*it)->production_vertex()->position().x()
		  << "\t" << (*it)->production_vertex()->position().y()
		  << "\t" << (*it)->production_vertex()->position().z()
		  << "\t" << (*it)->production_vertex()->position().t();
      std::cout << std::endl;
    }

    gen_n++;
  }
  
  return true;
}

int GlobeGenerator::mother(const HepMC::GenEvent* g, HepMC::GenParticle *p) {
  
  if (!(p->production_vertex())) 
    return 0;

  HepMC::GenVertex* inVertex = p->production_vertex();
  for(HepMC::GenVertex::particles_in_const_iterator iter = inVertex->particles_in_const_begin();
      iter != inVertex->particles_in_const_end(); ++iter) {
    
    int index = 0;
    for (HepMC::GenEvent::particle_const_iterator it = g->particles_begin(); it != g->particles_end(); ++it) {
      if (*it == *iter)
        return index;
      index++;
    }
  }
  
  return -1;
}
