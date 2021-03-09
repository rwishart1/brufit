#include "WignerGen.h"
//#include "lsconstants.h"

namespace HS{
  namespace FIT{

void WignerGen::do_line0 (double *l1, int j)
   {
   double xj = 1./j;
   l1[j] = -p*l1[j-1];
   for (int i=j-1; i>=1; --i)
     l1[i] = xj*sqt[j]*(q*sqt[j-i]*l1[i] - p*sqt[i]*l1[i-1]);
   l1[0] = q*l1[0];
   }


 void WignerGen::do_line (const double *l1, double *l2, int j, int k)
   {
   double xj = 1./j;
   double t1 = xj*sqt[j-k]*q, t2 = xj*sqt[j-k]*p;
   double t3 = xj*sqt[k]*p,   t4 = xj*sqt[k]*q;
   l2[j] = sqt[j] * (t4*l1[j-1]-t2*l2[j-1]);
   for (int i=j-1; i>=1; --i)
     l2[i] = t1*sqt[j-i]*l2[i] - t2*sqt[i]*l2[i-1]
             +t3*sqt[j-i]*l1[i] + t4*sqt[i]*l1[i-1];
   l2[0] = sqt[j] * (t3*l1[0]+t1*l2[0]);
   }
 
    WignerGen::WignerGen(int lmax, double ang)/*:p(TMath::Sin(ang/2)), q(TMath::Cos(ang/2)), sqt(2*lmax+1),
				     d(lmax+1,2*lmax+1), n(-1)*/
    { for (tsize m=0; m<sqt.size(); ++m) sqt[m] = TMath::Sqrt(double(m)); }


 arr2<double> &WignerGen::recurse ()
   {
   ++n;
   if (n==0)
     d[0][0] = 1;
   else if (n==1)
     {
     d[0][0] = q*q; d[0][1] = -p*q*sqt[2]; d[0][2] = p*p;
     d[1][0] = -d[0][1]; d[1][1] = q*q-p*p; d[1][2] = d[0][1];
     }
   else
     {
     int sign = (n&1)? -1 : 1;
     for (int i=0; i<=2*n-2; ++i)
       {
       d[n][i] = sign*d[n-2][2*n-2-i];
       sign=-sign;
       }
     do_line (d[n-1],d[n],2*n-1,n);
     for (int k=n; k>=2; --k)
       {
       do_line (d[k-2],d[k-1],2*n-1,k-1);
       do_line (d[k-1],d[k],2*n,k);
       }
     do_line0 (d[0],2*n-1);
     do_line (d[0],d[1],2*n,1);
     do_line0 (d[0],2*n);
     }
   return d;
   }
 

  }//FIT
}//HS
