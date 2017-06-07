#ifndef TUNNEL_CC
#define TUNNEL_CC

#include "vector.h"
#include "fresnel.h"
#include "resutil.h"
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

complex<double> tunnel (Vector<double> P, GlobalParms parms,
              double n1, double n2) 
{
 double cosbeta,rci,k;
 complex<double> Erg;
 rci=real(abs(P)/parms.n0);
 cosbeta=rci/parms.r0;
 
 k=parms.r0*sqrt(n2*n2/(n1*n1)*cosbeta*cosbeta-1.0);
 Erg=exp(-k*(parms.r0-abs(P)));
 return Erg; 
}


complex<double> tunnel_reflect (int pol, GlobalParms parms, Vector<double> P, Vector<double> k, double n1, double n2)
{
 Vector<double> n;
 double alpha,gamma,cosa,absrF2,abstF2,b,d;
 complex<double>  beta,rF,tF,Erg,rK,x;

 n=P/abs(P);
 alpha=acos(n*k);
 if (alpha>M_PI/2.0) alpha=M_PI-alpha;
 cosa=cos(alpha);
 beta=asin((complex<double> )(n1/n2*sin(alpha)));
 gamma=real(I*parms.k0*cos(beta));
 tF=Fresnel_trans(pol,k,n,n1,n2);
 rF=Fresnel_reflect(pol,k,n,n1,n2);
 abstF2=real(tF*conj(tF));
 absrF2=real(rF*conj(rF));
  b=cos(M_PI/2.0-alpha)*parms.r0*n1/n2;
 d=b-parms.r0;
 x=cosa-abstF2*exp(2.0*gamma*d); 
 rK=sqrt((complex<double>)(cosa-abstF2*exp(2.0*gamma*d))/(complex<double>)(absrF2*cosa));
 Erg=rK;
 return Erg;
}
#endif 

