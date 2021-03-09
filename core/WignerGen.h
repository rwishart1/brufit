//#pragma once
#ifndef PLANCK_WIGNER_H
#define PLANCK_WIGNER_H

#include <cmath>
#include "TMath.h"
#include "arr.h"
#include "sse_utils_cxx.h"

namespace HS{
  namespace FIT{

    class WignerGen{

    public:
      /*
	WignerGen()=default;
	WignerGen(const WignerGen&)=default;
	WignerGen(WignerGen&&)=default;
	virtual ~WignerGen()=default;
	WignerGen& operator=(const WignerGen& other)=default;
	WignerGen& operator=(WignerGen&& other)=default; 
      */
      WignerGen(int lmax, double ang)/*:p(TMath::Sin(ang/2)), q(TMath::Cos(ang/2)), sqt(2*lmax+1),
				       d(lmax+1,2*lmax+1), n(-1){}*/;
      
      arr2<double> &recurse ();
      
    private:

      double p,q;
      arr<double> sqt;
      arr2<double> d;
      int n;
      
      void do_line0 (double *l1, int j);
      void do_line (const double *l1, double *l2, int j, int k);
      
      //ClassDefOverride(HS::FIT::WignerGen, 1);
    };
    
  }//FIT
}//HS

#endif
