
/// Based on AmpTool BreitWigner amplitude

#include "RelBreitWigner.h"
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <TMath.h>
#include <cmath>

namespace HS{
  namespace FIT{

    RelBreitWigner::RelBreitWigner(const char *name, const char *title,
				   RooAbsReal& _x,
				   RooAbsReal& _m1,
				   RooAbsReal& _m2,
				   RooAbsReal& _L,
				   RooAbsReal& _mean,
				   RooAbsReal& _width) :
      HS::FIT::RooHSEventsPDF(name,title),
      x("x","x",this,_x),
      m1("m1","m1",this,_m1),
      m2("m2","m2",this,_m2),
      L("L","L",this,_L),
      mean("mean","mean",this,_mean),
      width("width","width",this,_width)
    {
      MakeSets();
      x.SetName(_x.GetName());
      m1.SetName(_m1.GetName());
      m2.SetName(_m2.GetName());
      L.SetName(_L.GetName());
      mean.SetName(_mean.GetName());
      width.SetName(_width.GetName());
    }


    RelBreitWigner::RelBreitWigner(const RelBreitWigner& other, const char* name) :
      HS::FIT::RooHSEventsPDF(other,name),
      x("x",this,other.x),
      m1("m1",this,other.m1),
      m2("m2",this,other.m2),
      L("L",this,other.L),
      mean("mean",this,other.mean),
      width("width",this,other.width)
    {
      MakeSets();
      x.SetName(other.x.GetName());
      m1.SetName(other.m1.GetName());
      m2.SetName(other.m2.GetName());
      L.SetName(other.L.GetName());
      mean.SetName(other.mean.GetName());
      width.SetName(other.width.GetName());
      if(fEvTree) SetEvTree(fEvTree,fCut);//Needs fProxSet filled first
    }
    void RelBreitWigner::MakeSets(){
      fProxSet.push_back(&x);
      fProxSet.push_back(&m1);
      fProxSet.push_back(&m2);
      fProxSet.push_back(&L);
      fParSet.push_back(&mean);
      fParSet.push_back(&width);
      InitSets();
    }



    Double_t RelBreitWigner::evaluate() const
    {
      Double_t mass  = x;

      // assert positive breakup momenta
      Double_t q0 = abs( breakupMomentum(mean, m1, m2) );
      Double_t q  = abs( breakupMomentum(mass, m1, m2) );

      Double_t F0 = barrierFactor(q0, L);
      Double_t F  = barrierFactor(q,  L);

      Double_t w = width*(mean/mass)*(q/q0)*((F*F)/(F0*F0));
      //Double_t w = width;

      // this first factor just gets normalization right for BW's that have
      // no additional s-dependence from orbital L
      complex<Double_t> bwtop( sqrt( mean * width / 3.1416 ) , 0.0 );

      complex<Double_t> bwbottom( ( mean*mean - mass*mass ) , -1.0 * ( mean * w ) );

      return(abs( F * bwtop / bwbottom ) * abs( F * bwtop / bwbottom ) );
    }

    Double_t RelBreitWigner::evaluateMC(const vector<Float_t> *vars,const  vector<Int_t> *cats) const {
      // ENTER IDENTICAL EXPRESSION TO evaluate() IN TERMS OF MC VARIABLE ARGUMENTS HERE
      Double_t mcx=(*vars)[fTreeEntry*fNvars+0];
      Double_t mcm1=(*vars)[fTreeEntry*fNvars+1];
      Double_t mcm2=(*vars)[fTreeEntry*fNvars+2];
      Double_t mcL=(*vars)[fTreeEntry*fNvars+3];
      Double_t mass  = mcx;

      // assert positive breakup momenta
      Double_t q0 = abs( breakupMomentum(mean, mcm1, mcm2) );
      Double_t q  = abs( breakupMomentum(mass, mcm1, mcm2) );

      Double_t F0 = barrierFactor(q0, mcL);
      Double_t F  = barrierFactor(q,  mcL);

      Double_t w = width*(mean/mass)*(q/q0)*((F*F)/(F0*F0));
      //Double_t w = width;

      // this first factor just gets normalization right for BW's that have
      // no additional s-dependence from orbital L
      complex<Double_t> bwtop( sqrt( mean * width / 3.1416 ) , 0.0 );

      complex<Double_t> bwbottom( ( mean*mean - mass*mass ) , -1.0 * ( mean * w ) );

      return(abs( F * bwtop / bwbottom ) * abs( F * bwtop / bwbottom ) );
    }


  }
}
