#ifndef FUNSURF_H
#define FUNSURF_H

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "funktion.h"
#include "form.h"
// #include "dreiecklist.h"


#ifndef SW
#define SW 20.0
#endif


/**
  *@author Thomas Weigel
  */

class Funsurf : public Form  {
public: 

   ~Funsurf();
   Funsurf();
   Funsurf(string &S);
   Funsurf (const Form &F);
   Funsurf (const Funsurf &FS);
   Funsurf ( const Funktion &Fkt,
             const Vector<double> &P,
             const Matrix <complex<double> > &alpha,
             const double r0,
             const Vector <double> &Ex,
             const Vector <double> &Ey,
             const Vector <double> &Ez,
             const complex<double>  n);

   Funsurf ( const string &S,
             const Vector<double> &P,
             const Matrix <complex<double> > &alpha,
             const double r0,
             const Vector <double> &Ex,
             const Vector <double> &Ey,
             const Vector <double> &Ez,
             const complex<double>  n);
	     
 bool next (const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside=-1);
 void scale(double sf){}
 double Volume () {return 0.0;}  // Muss noch angepasst werden !!!
 string fkt() {return Fkt.getfktm(); }
 bool  nullstelle(Funktion f,double &ns);
 bool isInside (const Vector<double> &p){return true;}
 Vector<double> norm (const Vector<double> &p) ;
 void setr0(double r0) {this->r0=r0;}
 void BerechneSurfGitter(int ntheta, int nphi);
 Vector<double> **SurfGitter;
 void initQuad (){};
 void setFkt (string f){Fkt=Funktion(f);}
 void binWrite (ofstream &os);
 void binRead(ifstream &os); 
 void setP (Vector<double> r) {P=r;initQuad();} // Setze Bezugspunkt des Einschlusses
  void setP (double x, double y, double z) {setP(Vector<double>(x,y,z)); }
 int nphi,ntheta;
protected :
  Funktion Fkt;
};

#endif

   
