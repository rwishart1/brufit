#pragma once

#include "RooHSComplex.h"
#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <Math/SpecFunc.h>
#include <TMath.h>
#include "arr.h"


namespace HS{
namespace FIT{

class RooHSDWigner : public RooHSComplex {

 public:
  RooHSDWigner()=default;
  RooHSDWigner(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, int l, int m=0,int s=0, double factor=1);
  RooHSDWigner(const RooHSDWigner& other, const char* name=nullptr);
  TObject* clone(const char* newname) const override{ return new RooHSDWigner(*this, newname);}
  ~RooHSDWigner() override=default;


TString FactoryReal() const final{
return TString("RooHSDWignerRe::")
+"Re"+GetName()+Form("(%s,%s,%s)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
}

TString FactoryImag() const final{
return TString("RooHSDWignerIm::")
+"Im"+GetName()+Form("(%s,%s,%s)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
}

TString FactoryImagConj() const final {
return TString("RooHSDWignerIm::")
+"CoIm"+GetName()+Form("(%s,%s,%s,-1)",_theta.arg().GetName(),_phi.arg().GetName(),GetName());
}


Int_t L() const{return _L;}
Int_t M() const{return _M;}
Int_t S() const{return _S;}

protected:

/*arr2<double>*/ void recurse () const;
Double_t evaluate() const final;

Bool_t CheckClean() const;
 
private:

RooRealProxy _theta;
RooRealProxy _phi;
mutable Double_t _p; //Sin(theta/2)
mutable Double_t _q; //Cos(theta/2)
Double_t _N=0;
mutable Double_t _lastTheta=1E10;
mutable Double_t _lastVal=0;
Int_t _L=0;
Int_t _M=0;
Int_t _S=0;
mutable Int_t _n;
mutable arr<double> _sqt;
mutable arr2<double> _dmatrix;

void do_line0 (double *l1, int j) const;
void do_line (const double *l1, double *l2, int j, int k) const;

      
ClassDefOverride(HS::FIT::RooHSDWigner,1);
};

  ////////////////////////////////////////////////////////////////////////////////
    
    inline Double_t RooHSDWigner::evaluate() const
    {
      if(CheckClean()) return _lastVal;
//CHANGE THIS TO DWIGNER EVALUATE!
//  return _lastVal=_N*ROOT::Math::assoc_legendre(_L,_absM,_lastCTheta);
return _lastVal = _dmatrix[_M][_S];
    }

////////////////////////////////////////////////////////////////////////////////

    inline Bool_t RooHSDWigner::CheckClean() const
    {
      Double_t th=_theta;
      if(th==_lastTheta ) return kTRUE;
      _lastTheta=th;
_p = TMath::Sin(th/2);
_q = TMath::Cos(th/2);
recurse();
      return kFALSE;
    }


///////////////////////////////////////////////////////////////////////////////////////////////
class RooHSDWignerIm : public RooAbsReal {
public:
RooHSDWignerIm() =default;
   
RooHSDWignerIm(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner, Int_t conj=1 );

RooHSDWignerIm(const RooHSDWignerIm& other, const char* name = nullptr);
TObject* clone(const char* newname) const override { return new RooHSDWignerIm(*this, newname); }
~RooHSDWignerIm() override =default;

protected:
      Double_t evaluate() const final;
      Bool_t CheckClean() const;

private:

      RooRealProxy _theta;
      RooRealProxy _phi;
      RooRealProxy _mag;
      mutable Double_t _lastTheta=-1E10;
      mutable Double_t _lastPhi=-1E10;
      mutable Double_t _lastVal=0;
      Int_t _M=0;
      Short_t _conj=1;

 ClassDefOverride(HS::FIT::RooHSDWignerIm,1);
    };

 ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerIm::evaluate() const{
      if(!_M)return 0; //M=0 =>real
      if(CheckClean()) return _lastVal;
//////CHANGE THIS FOR WIGNER FUNC
      Double_t angle=_M*_lastPhi; 
      return _lastVal=_conj*_mag*TMath::Sin(angle);
    }

   ////////////////////////////////////////////////////////////////////////////////
    
    inline Bool_t RooHSDWignerIm::CheckClean() const
    {
      Double_t th=_theta;
      Double_t ph=_phi;
      if(th==_lastTheta&&ph==_lastPhi ) return kTRUE;
      _lastTheta=th;
      _lastPhi=ph;
return kFALSE;
    }


 ///////////////////////////////////////////////////////////////////////////////////////////////
    class RooHSDWignerRe : public RooAbsReal {
    public:
      RooHSDWignerRe() =default;
      RooHSDWignerRe(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner );

      RooHSDWignerRe(const RooHSDWignerRe& other, const char* name = nullptr);
      TObject* clone(const char* newname) const override { return new RooHSDWignerRe(*this, newname); }
      ~RooHSDWignerRe() override =default;

 protected:
      Double_t evaluate() const final;
      Bool_t CheckClean() const;
    private:

RooRealProxy _theta;
      RooRealProxy _phi;
      RooRealProxy _mag;
      mutable Double_t _lastTheta=-1E10;
      mutable Double_t _lastPhi=-1E10;
      mutable Double_t _lastVal=0;
      Int_t _M=0;
 
      ClassDefOverride(HS::FIT::RooHSDWignerRe,1);
    };

 ////////////////////////////////////////////////////////////////////////////////
    inline Double_t RooHSDWignerRe::evaluate() const{
      if(CheckClean()) return _lastVal;
////CHANGE FOR WIGNER FUNC
      Double_t angle=_M*_lastPhi;
      return _lastVal=_mag*TMath::Cos(angle);
    }

////////////////////////////////////////////////////////////////////////////////
    
    inline Bool_t RooHSDWignerRe::CheckClean() const
    {
      Double_t th=_theta;
      Double_t ph=_phi;
      if(th==_lastTheta&&ph==_lastPhi ) return kTRUE;
      _lastTheta=th;
      _lastPhi=ph;
      return kFALSE; 
    }

}//FIT
}//HS
