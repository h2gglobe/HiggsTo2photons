#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
 
class MyEventFilter : public edm::EDFilter {
public:
  MyEventFilter(const edm::ParameterSet&);
  ~MyEventFilter();
private:
  bool filter( edm::Event &, edm::EventSetup const& );
  std::vector<edm::EventRange> vEventRange_;
  
};
 
#include <iostream>
 
using namespace std;
using namespace edm;
 
MyEventFilter::MyEventFilter(const ParameterSet & pset) :
  vEventRange_(pset.getUntrackedParameter<std::vector<edm::EventRange> >("vEventRange")) 
{}
 
MyEventFilter::~MyEventFilter() 
{}
 
bool MyEventFilter::filter(edm::Event &e, edm::EventSetup const& es) {

  //std::cout << e.id().event() << " " << e.id().run() << std::endl;
  for (unsigned int i=0; i<vEventRange_.size(); i++) {
    edm::EventRange range = vEventRange_[i];

    if((e.id().event() >= range.startEvent() and e.id().event() <= range.endEvent()) and
       (e.id().run() >= range.startRun() and e.id().run() <= range.endRun())) {
      //std::cout << "PASS" << std::endl;
      return true;
    }
  }

  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MyEventFilter);
