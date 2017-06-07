/***************************************************************************
                          istrahl.h  -  description                              
                             -------------------                                         
    begin                : Fri Oct 15 1999                                           
    copyright            : (C) 1999 by Thomas Weigel                         
    email                : weigel@lat.ruhr-uni-bochum.de                                     
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   * 
 *                                                                         *
 ***************************************************************************/


#ifndef ISTRAHL_H
#define ISTRAHL_H
#include <complex>
#include "fresnel.h"
#include "resutil.h"
#include "vector.h"
#include "matrix.h"
#include "ebene.h"
#include "form.h"
#include "strahl.h"
#include "misc.h"
// #define KORR 1.0e-8



/**
  Klasse zur Strahlberechnung, mit zwei elektrischen Feldern, deren Polarisationsrichtung senkrecht zueinander stehen. Diese Klasse dient zur simultanen Berechnung der E-Felder für beide Polarisationsrichtungen
  *@author Thomas Weigel
  */
 

class IStrahl {

public: 
  IStrahl();
  IStrahl(const Vector<double> &p,
         const Vector<complex<double> > &Pol, const Vector<double> &K,
         complex<double>  n0, double r0, double k0,
         const int Anzein=0, Form **Einschluss=NULL);
  IStrahl(const IStrahl &r)
  {
	  this->E1=r.E1;
	  this->E2=r.E2;
	 
	  this->AnzEin=r.AnzEin;
	  copyFormList(this->Ein,r.Ein,r.AnzEin);
	  // this->Ein=r.Ein; 
	  this->einindex=r.einindex;
	  this->getunnelt=r.getunnelt;
	  this->imEinschluss=r.imEinschluss;
	  this->iR=r.iR;
	  this->isValid=r.isValid;
	  this->k=r.k;
	  this->k0=r.k0;
	  this->KORR=r.KORR;
	  this->logRay=r.logRay;
	  this->n=r.n;
	  this->OK=r.OK;
	  this->P=r.P;
	  this->r0=r.r0;	  
  }


  void next();
  IStrahl reflect(Vector<double> n, complex<double>  n1, complex<double>  n2);
  void refract(Vector<double> N, complex<double>  n1, complex<double>  n2);
  void tunnel(Vector<complex<double> > Pol, complex<double>  n1, complex<double>  n2);
  void tunnel(Vector<complex<double> > Pol, complex<double>  n1, complex<double>  n2, int l);
  double flaeche();
  void setRefract (complex<double>  n) { this->n=n; }  
  complex<double>  getRefract() { return n; }
  void setk(const Vector<double> &K) { k=K; }
  void setP(const Vector<double> &p) { P=p; }
  void setiR (int i) { iR=i; }
  void setGetunnelt (bool v) { getunnelt=v; }
  Vector<double> getk() { return k; }
  Vector<double> getP() { return P; }


  Form* Einschluss(int i) {return Ein[i]; }
  ~IStrahl();
  bool Einschluss() {return imEinschluss;}
  void checkEinschluss(int &Index, Vector<double> &Pmin);
  // Vector<double> checkEinschluss(const Vector<double>& Ps,int& Index);
  bool Getunnelt() { return getunnelt; }
  int EinIndex() { return einindex; }

  friend ostream& operator << (ostream &os,IStrahl S);
  int reflections () {return iR;}


public:
	Vector<complex<double> > Pol1() {return E1/abs(E1); }
	Vector<complex<double> > Pol2() {return E2/abs(E2); }
  void setKorr (double Korr) {KORR=Korr; }
  double  getKorr () {return KORR; } 
  Matrix<complex<double> > Fresnel_reflect (double alpha,complex<double>  n1, complex<double>  n2);
  Matrix<complex<double> > Fresnel_trans (double alpha,complex<double>  beta, complex<double>  n1, complex<double>  n2);
  void init_Efeld (const Ebene& Eb, const Vector<complex<double> >& Pol,const int AnzRays=1);
  void init_Efeld (const Vector<complex<double> >& PolS, const Vector<complex<double> > &PolP,
                         const int AnzRays);
  void init_Efeld (const Ebene& Eb, const Vector<complex<double> >& Pol1, const Vector<complex<double> >& Pol2, const int AnzRays);
  void initEfeldFokus (double sigma2, Vector<double> focuspos,  Vector<complex<double> > Pol);
void init_EfeldGauss (const Ebene& Eb,
                                   const Vector<complex<double> >& PolS,
                                   const Vector<complex<double> > &PolP,
                                   Gauss g);

  void initGauss (Vector<complex<double> > &Pol, Gauss g);
  double cross (const Vector<double> P10, const Vector<double> P11,
                       const Vector<double> P20, const Vector<double> P21);
  double crossXAxis (const Vector<double>& P1, const Vector<double>& P2,const Vector<double> &k);
  Vector<double> crossPlane (const Vector<double> Pe, const Vector<double> n);
  Vector<double> intersectRect (const Vector<double> P, const Vector<double> e1, const Vector<double> e2); /// Berechnung der Schnittpunkte mit einem Rechteck P:Ecke, e1,e2 : Vektoren entlang zweier Kanten die von P loslaufen. Gibt einen NaN-Vektor zurück, falls es keinen Schnittpunkt gibt

 // project(int i){};
  Vector<double> OK;
  Vector<double> P;
  Vector<double> k;
  Vector<complex<double> > E1,E2;
  Form **Ein;
  int AnzEin;
  int einindex;
  int iR;
  bool imEinschluss;
  bool getunnelt;
  bool logRay;
  complex<double>  n,k0;
  double KORR;
  double r0;
  bool isValid;
};


#endif

