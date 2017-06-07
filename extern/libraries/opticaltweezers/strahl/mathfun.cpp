// Routinen zur Berechnung v. Bessel-, Neumann- und Hankelfunktionen
// (hoffentlich) durch Rückwärtsrekursion
#include<complex>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>

using namespace std;

double envj(int n, double x)
{
 double env = 0.5*log10(6.28*n)- n*log10(1.36*x/n);
 return env;
}

int msta1(double x, int mp)
{
 double a0, f0, f1, f;
 int n0, n1, nn;
 a0 = fabs(x);
 n0 = int(1.1*a0)+1;
 f0 = envj(n0,a0)-mp;
 n1 = n0+5;
 f1 = envj(n1,a0)-mp;
 for(int it=1;it<=20;it++)
 {
  nn = int(n1-(n1-n0)/(1.0-f0/f1));
  f = envj(nn,a0)-mp;
  if(fabs((double)(nn-n1))<1.0)
   break;
  n0 = n1;
  f0 = f1;
  n1 = nn;
  f1 = f;
 }
// cout << "nn:" << nn << "\n";
 return nn;
}

int msta2(double x, int n, int mp)
{
 double a0, hmp, ejn, f0, f1, f, obj;
 int n0, n1, nn;
 a0 = fabs(x);
 hmp = 0.5*mp;
 ejn=envj(n,a0);
 if(ejn<=hmp)
 {
  obj = mp;
  n0 = int(1.1*a0);
//  cout <<"n0" << n0 << "\n";
 }
 else
 {
  obj = hmp + ejn;
  n0 = n;
//  cout <<"n0" << n0 << "\n";
 }
 f0 = envj(n0,a0)-obj;
 n1 = n0+5;
 f1 = envj(n1,a0)-obj;
 for(int it=1;it<=20;it++)
 {
  nn = n1 - (n1 - n0)/(1.0-f0/f1);
//  cout << "nn1:" << nn << "\n";
  f = envj(nn,a0)-obj;
  if (abs(nn-n1)<1)
    break;
  n0 = n1;
  f0 = f1;
  n1 = nn;
  f1 = f;
 }
// cout << "nn+10:" << nn+10 << "\n";
 return (nn+10);
}

int jyndif(int nmax, complex<double>  z, complex<double>  *j, complex<double>  *jd, 
           complex<double>  *y, complex<double>  *yd)
{
 double a0; 
 complex<double>  csa, csb, cs, cf0, cf1, cf, *jh;
 int m, ee;

 a0 = abs(z);

// argument == 0 ?
 if(a0<1.0e-60)
 {
  for (int i=0;i<=nmax;i++)
  {
    j[i] = complex<double> (0.0,0.0);
   jd[i] = complex<double> (0.0,0.0);
    y[i] = complex<double> (-1.e300,0.0);
   yd[i] = complex<double> (1.0e300,0.0);
  }
   j[0] = complex<double> (1.0,0.0);
   jd[1] = complex<double> (1.0/3.0,0.0);
  return 0;
 }

 j[0] = sin(z)/z;
 j[1] = (j[0]-cos(z))/z;

 if(nmax>=2)
 {
   csa=j[0];
   csb=j[1];
   m = msta1(a0,200);
   if (m<nmax)
   {
     nmax = m;
     ee = 1;
   }
   else
   {
     m =  msta2(a0,nmax,17);
     ee = 0;
   }
   jh = new complex<double> [m+1];
   cf0 = complex<double> (0.0,0.0);
   cf1 = complex<double> (1.0e0-100,0.0);

   for(int k=m;k>=0;k--)   
   {
     cf = (2.0*k+3.0)*cf1/z - cf0;
     if (k<=nmax)
       jh[k]=cf;
     cf0 = cf1;
     cf1 = cf;
   }
   if((abs(csa))>(abs(csb)))
      cs = csa / cf;
   if((abs(csa))<=(abs(csb)))
      cs = csb / cf0;
   for(int k=0;k<=nmax;k++)
   {
     j[k] = cs*jh[k];   
   }
   delete[] jh;
 } 
 jd[0] = (cos(z)-sin(z)/z)/z;
 for(int k=1;k<=nmax;k++)
 {
  jd[k] = j[k-1]-(k+1.0)*j[k]/z; 
 }
 
 y[0] = -cos(z)/z;
 y[1] = (y[0]-sin(z))/z;
 yd[0] = (sin(z)+cos(z)/z)/z;
 yd[1] = (2.0*yd[0]-cos(z))/z;

 for(int k=2;k<=nmax;k++)
 {
   if(abs(j[k-1])>abs(j[k-2]))
     y[k] = (j[k]*y[k-1]-1.0/(z*z))/j[k-1];
   else
     y[k] = (j[k]*y[k-2]-(2.0*k-1)/(z*z*z))/j[k-2];
 }
 for(int k=2;k<=nmax;k++)
 {
   yd[k] = y[k-1] - (k+1.0)*y[k]/z;
 }
 return ee;
}


int hndif(int nmax, int type, complex<double>  z, complex<double>  *h, 
          complex<double>   *hd)
{
 int ee;
 complex<double>  *j, *jd, *y, *yd, ii;
 ii = complex<double> (0.0,1.0);
 j = new complex<double> [nmax+1];
 jd = new complex<double> [nmax+1];
 y = new complex<double> [nmax+1];
 yd = new complex<double> [nmax+1];
 ee=jyndif(nmax, z, j, jd, y, yd);

 for(int k=0;k<=nmax;k++)
 {
   if(type==1)
   {
     h[k]  = j[k]  + ii*y[k];
     hd[k] = jd[k] + ii*yd[k];
   }
   else 
     if(type==2)
     {
       h[k]  = j[k]  - ii*y[k];
       hd[k] = jd[k] - ii*yd[k];
     }
     else
     {
       cout << "Fehler, existiert nicht!\n";
       return 1;
     }
   }
 delete j;
 delete jd;
 delete y;
 delete yd;

 return 0; 
}

int jyndif(int nmax, double x, double *j, double *jd, double *y, double *yd)
{
 complex<double>  *jh, *jdh, *yh, *ydh;
 complex<double>  z = complex<double> (x,0.0);
 jh = new complex<double> [nmax+1];
 jdh = new complex<double> [nmax+1];
 yh = new complex<double> [nmax+1];
 ydh = new complex<double> [nmax+1];
 int ee =jyndif(nmax, z, jh, jdh, yh, ydh); 
 for(int k=0;k<=nmax;k++)
 {
   j[k]  = real(jh[k]);
   jd[k] = real(jdh[k]);
   y[k] = real(yh[k]);
   yd[k] = real(ydh[k]);
 }
 delete jh;
 delete jdh;
 delete yh;
 delete ydh;

 return ee;
}

int ricjyndif(int nmax, complex<double>  z, complex<double>  *j, complex<double>  *jd,
            complex<double>  *y, complex<double>  *yd)
{
 int ee = jyndif(nmax, z, j, jd, y, yd);

 for(int k=0;k<=nmax;k++)
 {
   j[k] = z*j[k];
   y[k] = z*y[k];
  jd[k] = jd[k]+j[k]/z;
  yd[k] = yd[k]+y[k]/z;
 }
 return ee;
}

int richndif(int nmax, int type, complex<double>  z, complex<double>  *h,
           complex<double>   *hd)
{
 int ee;
 complex<double>  *j, *jd, *y, *yd, ii;
 ii = complex<double> (0.0,1.0);
 j = new complex<double> [nmax+1];
 jd = new complex<double> [nmax+1];
 y = new complex<double> [nmax+1];
 yd = new complex<double> [nmax+1];
 ee=jyndif(nmax, z, j, jd, y, yd);

 for(int k=0;k<=nmax;k++)
 {
   if(type==1)
   {
     h[k]  = z*j[k]  + ii*z*y[k];
     hd[k] = jd[k] + j[k]/z + ii*(yd[k] + y[k]/z);
   }
   else
     if(type==2)
     {
       h[k]  = z*j[k]  - ii*z*y[k];
       hd[k] = jd[k] + j[k]/z - ii*(yd[k] + y[k]/z);
     }
     else
     {
       cout << "Fehler, existiert nicht!\n";
       return 1;
     }
   }
 delete j;
 delete jd;
 delete y;
 delete yd;

 return 0;
}

int ricjyndif(int nmax, double x, double *j, double *jd, double *y, double *yd)
{
 complex<double>  *jh, *jdh, *yh, *ydh;
 complex<double>  z = complex<double> (x,0.0);
 jh = new complex<double> [nmax+1];
 jdh = new complex<double> [nmax+1];
 yh = new complex<double> [nmax+1];
 ydh = new complex<double> [nmax+1];
 int ee = jyndif(nmax, z, jh, jdh, yh, ydh);

 for(int k=0;k<=nmax;k++)
 {
   j[k] = real(z*jh[k]);
   y[k] = real(z*yh[k]);
  jd[k] = real(jdh[k]+jh[k]/z);
  yd[k] = real(ydh[k]+yh[k]/z);
 }
 delete jh;
 delete jdh;
 delete yh;
 delete ydh;

 return ee;
}

