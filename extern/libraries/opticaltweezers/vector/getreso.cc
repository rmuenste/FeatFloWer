#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <complex>
#include "koef.h"
#include "getreso.h"

#define inclGH false
#define _USE_MATH_DEFINES 
#include <math.h>


#ifndef MAX_ANZ_VERSUCHE
#define MAX_ANZ_VERSUCHE 1000
#endif
using namespace std;

double F (int pol, int b, double n, double x, double l)
{
 double nx2, phi, ksi, delta;
 double eps, p,ls,ls2, shift, FH;
 double deltah;
 if (pol==TE) eps=1.0; 
 else eps=1/(n*n);
 p=b-1;
 ls=l+0.5; 
 //ls=sqrt(l*(l+1.0));
 ls2=ls*ls;
 nx2=n*n*x*x;
 
 ksi=sqrt(nx2-ls2);
 if (ls > x) 
 { 
  deltah=ls2-x*x;
  delta=sqrt(deltah);
    shift=atan(delta / eps / ksi);
 }
 else
  shift=0.0;
 phi=ksi - ls * acos(ls / n / x);
 
 FH=(p+0.25) * M_PI + shift - phi;
 if (inclGH)
  return FH - 0.5 * (sqrt(1.0-ls2/(n*n)/(x*x))-(ls / n / x)) * 10.0; 
 else
  return FH;
}

double Fneu(int pol, int b, double n, double x, double l)
{
 double nx2, phi, ksi, delta;
 double eps, p,ls,ls2, shift, FH;
 double deltah;
 complex<double>  rih;
 
 rih =  ri(l,  pol, complex<double> (1.0,0.0), complex<double> (n,0.0), x);
// cout << "rih:" << rih << "\n";
 if (pol==TE) eps=1.0;
 else eps=1/(n*n);
 p=b-1;
 ls=l+0.5;
 //ls=sqrt(l*(l+1.0));
 ls2=ls*ls;
 nx2=n*n*x*x;

 ksi=sqrt(nx2-ls2);
 if (ls > x)
 {
  deltah=ls2-x*x;
  delta=sqrt(deltah);
//    shift=atan(delta / eps / ksi);
//      shift = arg(rih)/2.0;
      shift=atan2(imag(rih),real(rih))/2.0;
//      cout << "shift:" << shift << "\n";
 }
 else
  shift=0.0;
 phi=ksi - ls * acos(ls / n / x);

 FH=(p+0.25) * M_PI + shift - phi;
 if (inclGH)
  return FH - 0.5 * (sqrt(1-ls2/(n*n)/(x*x))-(ls / n / x)) * 10;
 else
  return FH;
}



double GeoX (int pol, int b, double n, double x, double l)
{
 double xtrial, deltax, F1, F2, dfdx, dx;
 
 xtrial=l/(n+1.0) * 2.0;
 dx=1E-10;
 do
 {
  F1=F(pol,b,n,xtrial,l);
  F2=F(pol,b,n,xtrial+dx,l);
  dfdx=(F2-F1)/dx;
  deltax=F1/dfdx;
  xtrial-=deltax;
 }
 while (fabs(deltax) >= 1E-8);
 return xtrial;
}

double GeoXneu (int pol, int b, double n, double x, double l)
{
 double xtrial, deltax, F1, F2, dfdx, dx;
 int c=0;
 xtrial=l/(n+1.0) * 2.0;
 dx=1E-10;
 do
 {
  F1=Fneu(pol,b,n,xtrial,l);
  F2=Fneu(pol,b,n,xtrial+dx,l);
  dfdx=(F2-F1)/dx;
//  cout << "dfdx:" << dfdx << "\n";
  deltax=F1/dfdx;
  xtrial-=deltax;
  c++;
 }
 while ((fabs(deltax) >= 1E-10) && (c<MAX_ANZ_VERSUCHE));
 return xtrial;
}


double GeoXneuErr (int pol, int b, double n, double x, double l, double &deltax)
{
 double xtrial, F1, F2, dfdx, dx;
 int c=0;
 xtrial=l/(n+1.0) * 2.0;
 dx=1E-10;
 do
 {
  F1=Fneu(pol,b,n,xtrial,l);
  F2=Fneu(pol,b,n,xtrial+dx,l);
  dfdx=(F2-F1)/dx;
//  cout << "dfdx:" << dfdx << "\n";
  deltax=F1/dfdx;
  xtrial-=deltax;
  c++;
 }
 while ((fabs(deltax) >= 1E-10) && (c<MAX_ANZ_VERSUCHE));
 deltax=fabs(deltax);
 return xtrial;
}
