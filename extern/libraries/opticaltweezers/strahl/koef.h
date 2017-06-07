#include<complex>
using namespace std;

#ifndef KOEF_H
#define KOEF_H
complex<double>  ri(int l, int pol, complex<double>  na, complex<double>  np, double x);
complex<double>  tie(int l, int pol, complex<double>  na, complex<double>  np, double x);
complex<double>  tei(int l, int pol, complex<double>  na, complex<double>  np, double x);
complex<double>  rineu(int l, int pol, complex<double>  na, complex<double>  np, double x);
complex<double>  tineu(int l, int pol, complex<double>  na, complex<double>  np,
                  double x);
complex<double>  teneu(int l, int pol, complex<double>  na, complex<double>  np,
                  double x);
#endif
