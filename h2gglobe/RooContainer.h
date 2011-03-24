#ifndef ROOCONTAINER
#define ROOCONTAINER

// RooFit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"

#include "TRandom3.h"

#include <map>
#include <vector>
// For now only conatins expnentials and real vars and datasets.
// Abstract to addPdf soon

class RooContainer {

  public:
   RooContainer(){}
   ~RooContainer(){};

   double rand_n();
   
   void AddRealVar(const char*,float,float);
   void AddRealVar(const char*,float,float,float);
   void AddGenericPdf(const char*,const char*,
		      std::vector<const char*> &, double norm=10);
   void ComposePdf(const char* , const char *
			     ,std::vector<const char*> &);

   void CreateDataSet(const char*);
   void SetRealVar(const char * name, float x);
   
   void FitToData(const char*,const char*);
   void FitToData(const char*,const char *
	         ,double, double, double, double);

   
  private:
   std::map<const char*, RooRealVar> m_real_var_;
   std::map<const char*, RooExtendPdf> m_exp_;
   std::map<const char*, RooGenericPdf> m_gen_;
   std::map<const char*, RooAddPdf> m_pdf_;

   std::map<const char*, float> m_var_min_;
   std::map<const char*, float> m_var_max_;

   RooDataSet *data_;
   

};


#endif
