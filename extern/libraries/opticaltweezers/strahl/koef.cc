#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "mathfun.h"
#include "koef.h"
#include "vector.h"
using namespace std;

complex<double>  ri(int l, int pol, complex<double>  na, complex<double>  np,
                  double x)
{
 int ee;
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, ri;

// cout <<"l:" << l << " pol:" << pol << " na:" << na << " np:" << np << " x:" << x << "\n"; 
 
 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

 if(pol==1) // TE
   alphaq = 1.0;
 else
  if(pol==0) // TM
    alphaq = np / na;
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
     cout << "Pol =" << pol << "\n";
  }


 ee = hndif(l, 1, na*x, h1a, hd1a);
 ee = hndif(l, 1, np*x, h1p, hd1p);
 ee = hndif(l, 2, na*x, h2a, hd2a);
 ee = hndif(l, 2, np*x, h2p, hd2p);

 complex<double>  h1ad = (na*h1a[l]) + (na*na*x*hd1a[l]);
 complex<double>  h2ad = (na*h2a[l]) + (na*na*x*hd2a[l]);

 complex<double>  h1pd = (np*h1p[l]) + (np*np*x*hd1p[l]);
 complex<double>  h2pd = (np*h2p[l]) + (np*np*x*hd2p[l]);


 ri = - h1p[l]/h2p[l]*(na*x*h2a[l]*h2pd - alphaq*alphaq*np*x*h2p[l]*h2ad)/
        (na*x*h2a[l]*h1pd - alphaq*alphaq*np*x*h1p[l]*h2ad);
// ri = - (na*x*h2a[l]*h2pd - alphaq*alphaq*np*x*h2p[l]*h2ad)/
//        (na*x*h2a[l]*h1pd - alphaq*alphaq*np*x*h1p[l]*h2ad);


 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;

 return ri;

}

complex<double>  tie(int l, int pol, complex<double>  na, complex<double>  np, 
                   double x)
{
 int ee;
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, tie;
 complex<double>   ii = complex<double> (0.0,1.0);

 if(pol==1) // TE
   alphaq = 1.0;
 else
  if(pol==0) // TM
    alphaq = np / na;
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
  }
 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

// ee = richndif(l, 1, n1*x, h11, hd11);
// ee = richndif(l, 1, n2*x, h12, hd12);
// ee = richndif(l, 2, n1*x, h21, hd21);
// ee = richndif(l, 2, n2*x, h22, hd22);

ee = hndif(l, 1, na*x, h1a, hd1a);
ee = hndif(l, 1, np*x, h1p, hd1p);
ee = hndif(l, 2, na*x, h2a, hd2a);
ee = hndif(l, 2, np*x, h2p, hd2p);

complex<double>  h1ad = (na*h1a[l]) + (na*na*x*hd1a[l]);
complex<double>  h2ad = (na*h2a[l]) + (na*na*x*hd2a[l]);

complex<double>  h1pd = (np*h1p[l]) + (np*np*x*hd1p[l]);
complex<double>  h2pd = (np*h2p[l]) + (np*np*x*hd2p[l]);


 tie = h2a[l]/h2p[l]*2*alphaq*na*ii/
       (na*x*h2a[l]*h1pd-alphaq*alphaq*np*x*h1p[l]*h2ad);

 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;

 return tie;
}

complex<double>  tei(int l, int pol, complex<double>  na, complex<double>  np,
                   double x)
{
 int ee;
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, tei;
 complex<double>   ii = complex<double> (0.0,1.0);

 if(pol==1) // TE
   alphaq = 1.0;
 else
  if(pol==0) // TM
    alphaq = np / na;
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
  }
 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

// ee = richndif(l, 1, n1*x, h11, hd11);
// ee = richndif(l, 1, n2*x, h12, hd12);
// ee = richndif(l, 2, n1*x, h21, hd21);
// ee = richndif(l, 2, n2*x, h22, hd22);

ee = hndif(l, 1, na*x, h1a, hd1a);
ee = hndif(l, 1, np*x, h1p, hd1p);
ee = hndif(l, 2, na*x, h2a, hd2a);
ee = hndif(l, 2, np*x, h2p, hd2p);

complex<double>  h1ad = (na*h1a[l]) + (na*na*x*hd1a[l]);
complex<double>  h2ad = (na*h2a[l]) + (na*na*x*hd2a[l]);

complex<double>  h1pd = (np*h1p[l]) + (np*np*x*hd1p[l]);
complex<double>  h2pd = (np*h2p[l]) + (np*np*x*hd2p[l]);


 tei = h1p[l]/h1a[l]*2*alphaq*np*ii/
       (na*x*h2a[l]*h1pd-alphaq*alphaq*np*x*h1p[l]*h2ad);

 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;
            
 return tei;
}

complex<double>  rineu(int l, int pol, complex<double>  na, complex<double>  np,
                  double x)
{
 int ee;
 // Hankelfunktionen: h[d][1,2][a,p]
 // [1,2]: Typ ; [a,p] Aussen, Partikel
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, ri, rineu, Zij;

 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

 if(pol==1) // TE
 {
   alphaq = 1.0;
// Zij = Z_2 / Z_1 = n_1 / n_2 f. \mu_1 = \mu_2 = 1 
   Zij = na / np;
 }
 else
  if(pol==0) // TM
  {
    alphaq = np / na;
// Zij = Z_1 / Z_2 = n_2 / n_1 f. \mu_1 = \mu_2 = 1
    Zij = np / na;
  }
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
  }


 ee = hndif(l, 1, na*x, h1a, hd1a);
 ee = hndif(l, 1, np*x, h1p, hd1p);
 ee = hndif(l, 2, na*x, h2a, hd2a);
 ee = hndif(l, 2, np*x, h2p, hd2p);

// double Gamma = sqrt(l*(l+1));
double Gamma = (l+0.5);

 complex<double>  Ira = 1/(sqrt(na*x*sqrt(na*x*na*x-Gamma)));
 complex<double>  Irp = 1/(sqrt(np*x*sqrt(np*x*np*x-Gamma)));

// complex<double>  h1ad = (na*h1a[l]) + (na*na*x*hd1a[l]);
// complex<double>  h2ad = (na*h2a[l]) + (na*na*x*hd2a[l]);

 complex<double>  h1ad = 1.0/Ira * (h1a[l]/(na*x) + hd1a[l]) +
                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h1a[l];
 complex<double>  h2ad = 1.0/Ira * (h2a[l]/(na*x) + hd2a[l]) +
                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h2a[l];
 complex<double>  h1pd = 1.0/Irp * (h1p[l]/(np*x) + hd1p[l]) +
                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h1p[l];
 complex<double>  h2pd = 1.0/Irp * (h2p[l]/(np*x) + hd2p[l]) +
                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h2p[l];

// complex<double>  h1pd = (np*h1p[l]) + (np*np*x*hd1p[l]);
// complex<double>  h2pd = (np*h2p[l]) + (np*np*x*hd2p[l]);


// ri = - h1p[l]/h2p[l]*(na*x*h2a[l]*h2pd - alphaq*alphaq*np*x*h2p[l]*h2ad)/
//        (na*x*h2a[l]*h1pd - alphaq*alphaq*np*x*h1p[l]*h2ad);
 rineu =  - h2p[l]/h1p[l]*(h1a[l]/Ira*h1pd-Zij*h1p[l]/Irp*h1ad)/
          (h1a[l]/Ira*h2pd-Zij*h2p[l]/Irp*h1ad); // *Irp/Ira

 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;

 return rineu;
}

complex<double>  tineu(int l, int pol, complex<double>  na, complex<double>  np,
                  double x)
{
 int ee;
 // Hankelfunktionen: h[d][1,2][a,p]
 // [1,2]: Typ ; [a,p] Aussen, Partikel
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, ri, tineu, Zij;

 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

 if(pol==1) // TE
 {
   alphaq = 1.0;
// Zij = Z_1 / Z_2 = n_2 / n_1 f. \mu_1 = \mu_2 = 1
   Zij = np / na;
 }
 else
  if(pol==0) // TM
  {
    alphaq = na / np;
// Zij = Z_2 / Z_1 = n_1 / n_2 f. \mu_1 = \mu_2 = 1
    Zij = na / np;
  }
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
  }


 ee = hndif(l, 1, na*x, h1a, hd1a);
 ee = hndif(l, 1, np*x, h1p, hd1p);
 ee = hndif(l, 2, na*x, h2a, hd2a);
 ee = hndif(l, 2, np*x, h2p, hd2p);

// double Gamma = sqrt(l*(l+1));
double Gamma = (l+0.5);

// complex<double>  Ira = 1/(sqrt(na*x*sqrt(na*x*na*x-Gamma)));
// complex<double>  Irp = 1/(sqrt(np*x*sqrt(np*x*np*x-Gamma)));
 complex<double>  Ira=1.0;
 complex<double>  Irp=1.0;

 complex<double>  h1ad = h1a[l]/(na*x) + hd1a[l];
 complex<double>  h2ad = h2a[l]/(na*x) + hd2a[l];
 complex<double>  h1pd = h1p[l]/(np*x) + hd1p[l];
 complex<double>  h2pd = h2p[l]/(np*x) + hd2p[l];

// complex<double>  h1ad = 1.0/Ira * (h1a[l]/(na*x) + hd1a[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h1a[l];
// complex<double>  h2ad = 1.0/Ira * (h2a[l]/(na*x) + hd2a[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h2a[l];
// complex<double>  h1pd = 1.0/Irp * (h1p[l]/(np*x) + hd1p[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h1p[l];
// complex<double>  h2pd = 1.0/Irp * (h2p[l]/(np*x) + hd2p[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h2p[l];


// ri = - h1p[l]/h2p[l]*(na*x*h2a[l]*h2pd - alphaq*alphaq*np*x*h2p[l]*h2ad)/
//        (na*x*h2a[l]*h1pd - alphaq*alphaq*np*x*h1p[l]*h2ad);
 tineu =  (h2p[l]/Irp*h1pd-h1p[l]/Irp*h2pd)/
          (h2p[l]/Irp*h1ad-Zij*h1a[l]/Ira*h2pd)*alphaq*np/na*h1a[l]/h1p[l]*Irp/Ira; // *Irp/Ira

 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;

 return tineu;
}


complex<double>  teneu(int l, int pol, complex<double>  na, complex<double>  np,
                  double x)
{
 int ee;
 // Hankelfunktionen: h[d][1,2][a,p]
 // [1,2]: Typ ; [a,p] Aussen, Partikel
 complex<double>  *h1a, *hd1a, *h1p, *hd1p, *h2a, *hd2a, *h2p, *hd2p;
 complex<double>  alphaq, ri, teneu, Zij;

 h1a = new complex<double> [l+1];
 h1p = new complex<double> [l+1];
 hd1a = new complex<double> [l+1];
 hd1p = new complex<double> [l+1];
 h2a = new complex<double> [l+1];
 hd2a = new complex<double> [l+1];
 h2p = new complex<double> [l+1];
 hd2p = new complex<double> [l+1];

 if(pol==1) // TE
 {
   alphaq = 1.0;
// Zij = Z_1 / Z_2 = n_2 / n_1 f. \mu_1 = \mu_2 = 1
   Zij = np / na;
 }
 else
  if(pol==0) // TM
  {
    alphaq = na / np;
// Zij = Z_2 / Z_1 = n_1 / n_2 f. \mu_1 = \mu_2 = 1
    Zij = na / np;
  }
  else
  {
     cout << "Fehler";
     alphaq = 0.0;
  }


 ee = hndif(l, 1, na*x, h1a, hd1a);
 ee = hndif(l, 1, np*x, h1p, hd1p);
 ee = hndif(l, 2, na*x, h2a, hd2a);
 ee = hndif(l, 2, np*x, h2p, hd2p);

// double Gamma = sqrt(l*(l+1));
// double Gamma = (l+0.5);

// complex<double>  Ira = 1/(sqrt(na*x*sqrt(na*x*na*x-Gamma)));
// complex<double>  Irp = 1/(sqrt(np*x*sqrt(np*x*np*x-Gamma)));
 complex<double>  Ira=1.0;
 complex<double>  Irp=1.0;

 complex<double>  h1ad = h1a[l]/(na*x) + hd1a[l];
 complex<double>  h2ad = h2a[l]/(na*x) + hd2a[l];
 complex<double>  h1pd = h1p[l]/(np*x) + hd1p[l];
 complex<double>  h2pd = h2p[l]/(np*x) + hd2p[l];

// complex<double>  h1ad = 1.0/Ira * (h1a[l]/(na*x) + hd1a[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h1a[l];
// complex<double>  h2ad = 1.0/Ira * (h2a[l]/(na*x) + hd2a[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Ira*Ira*Ira * na*x*h2a[l];
// complex<double>  h1pd = 1.0/Irp * (h1p[l]/(np*x) + hd1p[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h1p[l];
// complex<double>  h2pd = 1.0/Irp * (h2p[l]/(np*x) + hd2p[l]) +
//                       (l+0.5)*(l+0.5)/2.0 *Irp*Irp*Irp * np*x*h2p[l];


// ri = - h1p[l]/h2p[l]*(na*x*h2a[l]*h2pd - alphaq*alphaq*np*x*h2p[l]*h2ad)/
//        (na*x*h2a[l]*h1pd - alphaq*alphaq*np*x*h1p[l]*h2ad);
 teneu =  (h1a[l]/Ira*h2ad-h2a[l]/Ira*h1ad)/
          (h1a[l]/Ira*h2pd-1.0/Zij*h2p[l]/Irp*h1ad)*na/np*1.0/alphaq*h2p[l]/h2a[l]*Ira/Irp; // *Irp/Ira

 delete h1a;
 delete h1p;
 delete hd1a;
 delete hd1p;
 delete h2a;
 delete hd2a;
 delete h2p;
 delete hd2p;

 return teneu;
}
