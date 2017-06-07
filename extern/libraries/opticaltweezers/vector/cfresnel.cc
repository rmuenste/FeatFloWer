#ifndef CFRESNEL_CC
#define CFRESNEL_CC

#include "fresnel.h" 
#include <fstream>
#include <complex>
#include <math.h>
 

complex<double>  Fresnel_trans (int Pol, complex<double>  alpha, double n1, double n2)

//  Berechnet den (komplexen) Fresnel-Koeffizienten fuer Transmission 
//  S : Strahl
//  n : Normale auf die reflektierende Stelle
//  n1 : Brechungsindex innen 
//  n2 : Brechungsindex aussen 
 {
   double n12;
   complex<double>  beta,Erg; 
  n12=n2/n1;
  beta=asin(sin(alpha) / n12);
  if (Pol==SENKRECHT) Erg=2.0*sin(beta)*cos(alpha)/sin(alpha+beta);
  else Erg=2.0 * sin(beta) * cos (alpha) / (sin(alpha+beta)*cos(alpha-beta));
  return Erg; 
}

#endif
