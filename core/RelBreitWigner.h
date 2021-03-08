#pragma once

#include <RooAbsPdf.h>
#include "RooHSEventsPDF.h"
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <complex>

namespace HS{
  namespace FIT{

    using std::complex;
    
    class RelBreitWigner : public RooHSEventsPDF {
  
    public:
      RelBreitWigner() = default; ; 
      RelBreitWigner(const char *name, const char *title,
		     RooAbsReal& _x,
		     RooAbsReal& _m1,
		     RooAbsReal& _m2,
		     RooAbsReal& _L,
		     RooAbsReal& _mean,
		     RooAbsReal& _width);
      RelBreitWigner(const RelBreitWigner& other, const char* name=nullptr) ;
      TObject* clone(const char* newname) const override { return new RelBreitWigner(*this,newname); }
      inline ~RelBreitWigner() override = default;

    protected:

      RooRealProxy x ;
      RooRealProxy m1 ;
      RooRealProxy m2 ;
      RooRealProxy L ;
      RooRealProxy mean ;
      RooRealProxy width ;
  
      Double_t evaluate() const override ;
      Double_t evaluateMC(const vector<Float_t> *vars,const  vector<Int_t> *cats) const override ;
      void MakeSets();

    private:
	
      // mass0 = mass of parent
      // mass1 = mass of first daughter
      // mass2 = mass of second daughter
      double breakupMomentum( double mass0, double mass1, double mass2 ) const {
	double q;
	// fabs -- correct?  consistent w/ previous E852 code
	q = sqrt( fabs(   mass0*mass0*mass0*mass0 + 
			  mass1*mass1*mass1*mass1 +
			  mass2*mass2*mass2*mass2 -
			  2.0*mass0*mass0*mass1*mass1 -
			  2.0*mass0*mass0*mass2*mass2 -
			  2.0*mass1*mass1*mass2*mass2  ) ) / (2.0 * mass0);
	return q;
      }

      // mass0 = mass of parent
      // spin  = angular momentum of the decay
      // mass1 = mass of first daughter
      // mass2 = mass of second daughter
      double barrierFactor( double mass0, int spin, double mass1, double mass2 ) const {
	double q;
	q = breakupMomentum(mass0, mass1, mass2);
	return barrierFactor( q, spin );
      }


      // q     = breakup momentum
      // spin  = angular momentum of the decay
      double barrierFactor ( double q, int spin ) const {
	double barrier;
	double z;
	z = ( (q*q) / (0.1973*0.1973) );
	switch (spin){
	case 0:
	  barrier = 1.0;
	  break;
	case 1:
	  barrier = sqrt( (2.0*z) / (z + 1.0) );
	  break;
	case 2:
	  barrier = sqrt( (13.0*z*z) / ((z-3.0)*(z-3.0) + 9.0*z) );
	  break;
	case 3:
	  barrier = sqrt( (277.0*z*z*z) / (z*(z-15.0)*(z-15.0)+9.0*(2.0*z-5.0)*(2.0*z-5.0)) );
	  break;
	case 4:
	  barrier = sqrt( (12746.0*z*z*z*z) / ((z*z-45.0*z+105.0)*(z*z-45.0*z+105.0)+25.0*z*(2.0*z-21.0)*(2.0*z-21.0)) );
	  break;
	default:
	  barrier = 0.0;
	}
	return barrier;
      }

      ClassDefOverride(HS::FIT::RelBreitWigner,1); // Your description goes here...
    };
  }
}
