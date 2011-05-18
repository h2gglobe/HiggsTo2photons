#include "SampleContainer.h"
#include <utility>
#include <iostream>

SampleContainer::~SampleContainer() 
{}

SampleContainer::SampleContainer() {
  weight = 1;
  itype = 0;
  ind = 0;
  histoplotit = 1;
  filesshortnam = "";
  ntot = 0;
  nred = 0;
  lumi = 0;
  xsec = 0;
  kfactor= 1;
  scale = 1;
  lumireal = 1;
  hasLumiSelection = false;
}

void SampleContainer::computeWeight(float intL) {
  if(itype==0) { //this is data
    weight = 1; 
  } else {
    weight = kfactor*scale*xsec*intL/ntot;
  }
}
 
// ----------------------------------------------------------------------------------------------------------------------
void SampleContainer::addGoodLumi(int run, int lumi1, int lumi2 )
{
	hasLumiSelection = true;
	goodLumis[run].push_back( std::make_pair(lumi1,lumi2) );
}

