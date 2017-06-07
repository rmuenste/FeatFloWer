#ifndef MATHFUN_H
#define MATHFUN_H 

#include "koef.h"

int jyndif(int nmax, complex<double>  z, complex<double>  *j, complex<double>  *jd,
            complex<double>  *y, complex<double>  *yd);
int hndif(int nmax, int type, complex<double>  z, complex<double>  *h, 
           complex<double>   *hd);
int jyndif(int nmax, double x, double *j, double *jd, double *y, double *yd);
int ricjyndif(int nmax, complex<double>  z, complex<double>  *j, complex<double>  *jd,
            complex<double>  *y, complex<double>  *yd);
int ricjyndif(int nmax, double x, double *j, double *jd, double *y, double *yd);
int richndif(int nmax, int type, complex<double>  z, complex<double>  *h,
           complex<double>   *hd);

#endif
