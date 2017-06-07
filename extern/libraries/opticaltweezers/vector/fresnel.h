#ifndef FRESNEL_H
#define FRESNEL_H
#include <complex>

#include "vector.h"
#include "resutil.h"


double abs (complex<double>  x);
double sqr (double x);
//dc Fresnel (StrahlInfo &S, Vector<double> n, double n1, double n2);
complex<double> Fresnel_trans (int pol, Vector<double> k, Vector<double> n, complex<double>  n1, complex<double>  n2);
complex<double> Fresnel_reflect (int pol, Vector<double> k, Vector<double> n, complex<double>  n1, complex<double>  n2);
complex<double>  freflect (int pol, double alpha, complex<double>  n1, complex<double>  n2);
complex<double>  ftrans (int pol, double alpha, complex<double>  n1, complex<double>  n2);
#endif
