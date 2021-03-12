#include "RooHSDWigner.h"

#include "TError.h"

namespace HS{
  namespace FIT{

void RooHSDWigner::do_line0 (double *l1, int j) const
   {
   double xj = 1./j;
   l1[j] = -_p*l1[j-1];
   for (int i=j-1; i>=1; --i)
     l1[i] = xj*_sqt[j]*(_q*_sqt[j-i]*l1[i] - _p*_sqt[i]*l1[i-1]);
   l1[0] = _q*l1[0];
   }


 void RooHSDWigner::do_line (const double *l1, double *l2, int j, int k) const
   {
   double xj = 1./j;
   double t1 = xj*_sqt[j-k]*_q, t2 = xj*_sqt[j-k]*_p;
   double t3 = xj*_sqt[k]*_p,   t4 = xj*_sqt[k]*_q;
   l2[j] = _sqt[j] * (t4*l1[j-1]-t2*l2[j-1]);
   for (int i=j-1; i>=1; --i)
     l2[i] = t1*_sqt[j-i]*l2[i] - t2*_sqt[i]*l2[i-1]
             +t3*_sqt[j-i]*l1[i] + t4*_sqt[i]*l1[i-1];
   l2[0] = _sqt[j] * (t3*l1[0]+t1*l2[0]);
   }

    RooHSDWigner::RooHSDWigner(const char* name, const char* title, RooAbsReal& theta, RooAbsReal& phi, int l, int m, int s,Double_t factor)
      : RooHSComplex(name, title)
      ,_theta("theta", "theta", this,theta)
      ,_phi("phi", "phi", this,phi)
      ,_sqt(2*l+1)
      ,_dmatrix(l+1,2*l+1) 
      ,_L(l),_M(m),_S(s)
      ,_n(-1)
    {
      for (tsize m=0; m<_sqt.size(); ++m) _sqt[m] = TMath::Sqrt(double(m));
    }

    /////////////////////////////////////////////////////////////////////////

    /*arr2<double>*/ void RooHSDWigner::recurse() const
   {
   ++_n;
   if (_n==0)
     _dmatrix[0][0] = 1;
   else if (_n==1)
     {
     _dmatrix[0][0] = _q*_q; 
     _dmatrix[0][1] = -_p*_q*_sqt[2]; 
     _dmatrix[0][2] = _p*_p;
     _dmatrix[1][0] = -_dmatrix[0][1]; 
     _dmatrix[1][1] = _q*_q-_p*_p; 
     _dmatrix[1][2] = _dmatrix[0][1];
     }
   else
     {
     int sign = (_n&1)? -1 : 1;
     for (int i=0; i<=2*_n-2; ++i)
       {
       _dmatrix[_n][i] = sign*_dmatrix[_n-2][2*_n-2-i];
       sign=-sign;
       }
     do_line (_dmatrix[_n-1],_dmatrix[_n],2*_n-1,_n);
     for (int k=_n; k>=2; --k)
       {
       do_line (_dmatrix[k-2],_dmatrix[k-1],2*_n-1,k-1);
       do_line (_dmatrix[k-1],_dmatrix[k],2*_n,k);
       }
     do_line0 (_dmatrix[0],2*_n-1);
     do_line (_dmatrix[0],_dmatrix[1],2*_n,1);
     do_line0 (_dmatrix[0],2*_n);
     }
   //  return d;
   }


  ////////////////////////////////////////////////////////////////////////////////

    RooHSDWigner::RooHSDWigner(const RooHSDWigner& other, const char* name)
      : RooHSComplex(other, name),
	_theta("theta", this, other._theta),
	_phi("phi", this, other._phi),
	//	_absM(other._absM),
	_L(other._L), _M(other._M)
    {
      // if(_M%2==1)_MFactor=-1;
      //else _MFactor =1;
      //_MFactor =1;
      _N=other._N;
    }


  ////////////////////////////////////////////////////////////////////////////////

    RooHSDWignerRe::RooHSDWignerRe(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner )
      :  RooAbsReal(name, title),
	 _theta("Theta","Theta",this,theta),
	 _phi("Phi","Phi",this,phi),
	 _mag("Mag","Mag",this,dwigner)
    {
      auto mag=dynamic_cast<RooHSDWigner*>(&dwigner);
      _M=mag->M();
    }

  ////////////////////////////////////////////////////////////////////////////////
    RooHSDWignerRe::RooHSDWignerRe(const RooHSDWignerRe& other, const char* name)
      : RooAbsReal(other, name),
	_theta("Theta", this, other._theta),
	_phi("Phi", this, other._phi),
	_mag("Mag", this, other._mag),
	_M(other._M)
    {
   
    }


 ////////////////////////////////////////////////////////////////////////////////

    RooHSDWignerIm::RooHSDWignerIm(const char *name, const char *title, RooAbsReal& theta, RooAbsReal& phi, RooAbsReal& dwigner,Int_t conj)
      :  RooAbsReal(name, title),
	 _theta("Theta","Theta",this,theta),
	 _phi("Phi","Phi",this,phi),
	 _mag("Mag","Mag",this,dwigner),
	 _conj(conj)
    {
      auto mag=dynamic_cast<RooHSDWigner*>(&dwigner);
      _M=mag->M();
    }

  ////////////////////////////////////////////////////////////////////////////////
    RooHSDWignerIm::RooHSDWignerIm(const RooHSDWignerIm& other, const char* name)
      : RooAbsReal(other, name),
	_theta("Theta", this, other._theta),
	_phi("Phi", this, other._phi),
	_mag("Mag", this, other._mag),
	_M(other._M),
	_conj(other._conj)
    {
   
    }
  
  }//FIT
}//HS
