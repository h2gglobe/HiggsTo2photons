#include "RooContainer.h"
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooFitResult.h>

#include <cmath>

using namespace RooFit;

void RooContainer::AddRealVar(const char* name ,float xmin,float xmax){
  RooRealVar temp(name,name,xmin,xmax);

  m_real_var_.insert(pair<const char*,RooRealVar>(name,temp));
  m_var_min_.insert(pair<const char*, float>(name,xmin));
  m_var_max_.insert(pair<const char*, float>(name,xmax));
  
  std::cout << "Added Variale"  << name << std::endl;
  
}

void RooContainer::AddRealVar(const char* name ,float init, float vmin, float vmax){
  RooRealVar temp(name,name,init,vmin,vmax);
  m_real_var_.insert(pair<const char*,RooRealVar>(name,temp));
  
}

void RooContainer::AddGenericPdf(const char* name,const char* formula, std::vector<const char*> & var){

    RooArgList roo_args;

    for (std::vector<const char*>::iterator it_var = var.begin()
	;it_var != var.end()
	;it_var++
	){
	  std::cout << "Adding Variable " << *it_var << std::endl;
	  roo_args.add(m_real_var_[*it_var]);
	}

    RooGenericPdf temp(name,name,formula,roo_args);
  			       
    m_exp_.insert(pair<const char*,RooGenericPdf>(name,temp));
}

void RooContainer::CreateDataSet(const char *name){
 
    data_ = new RooDataSet("data","data",RooArgSet(m_real_var_[name]) );

}

void RooContainer::FitToData(const char* name_func, const char * name_var){
 
    RooFitResult *fit_result = m_exp_[name_func].fitTo(*data_);

    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];
 
    RooPlot *xframe = m_real_var_[name_var].frame(x_min,x_max);
    data_->plotOn(xframe);
    m_exp_[name_func].plotOn(xframe);
 
    xframe->Draw();
}

void RooContainer::SetRealVar(const char * name, float x){

  
  std::map<const char*, RooRealVar>::const_iterator it_var  = m_real_var_.find(name);

    float min_x = m_var_min_[name];
    float max_x = m_var_max_[name];

    if (x > min_x && x < max_x){
      m_real_var_[name] = x;
      data_->add(RooArgSet(m_real_var_[name]));
    }

}

