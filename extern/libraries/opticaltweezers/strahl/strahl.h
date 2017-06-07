/***************************************************************************
                          strahl.h  -  description                              
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


#ifndef STRAHL_H
#define STRAHL_H
#include <complex>
//#include "ausgabefenster.h"
#include "fresnel.h"
#include "resutil.h"
#include "vector.h"
#include "matrix.h"
#include "ebene.h"
#include "form.h"


/*#ifndef KORR
#define KORR 1E-8
#endif */

#define KEIN_SCHNITTPUNKT     -1
#define GESCHNITTEN            0
#define IN_EBENE               1

#ifndef EPS_WINKEL 
#define EPS_WINKEL 0 /// Bei Brechung werden nur Strahlen weiterverfolgt, die einen Winkel > EPS_WINKEL gegenüber der Oberfläche aufweisen
#endif 
extern GlobalParms Parms;
extern double Grenze;
//extern AusgabeFenster *ausgabe;

/**
  *@author Thomas Weigel
  */


typedef struct StrahlBuffer
{
 /**
   Uebergabepuffer fuer den Konstruktor der Klasse >Strahl< zur Vereinfachung
 */
 complex<double>  n;
 Vector<double> P[5];
 Vector<double> k[5];
 Vector<complex<double> > E[5];
 int AnzEin;
 Form **Ein;
 double r0;
 Vector<complex<double> > ***Gitter;
};


/*#ifndef COMPLEX_I
#define COMPLEX_I
const complex<double>  I=complex<double>  (0.0,1.0);
#endif
  */

/// Parameter fuer Gaussstrahl
typedef struct
{
 double w0;          /// Strahlbreite im Fokus
 Vector<double> F;   /// Ort des Fokus
 Vector<double> k;   /// Strahlrichtung
 bool isGauss;
} Gauss;


void binWrite (Gauss gs, ofstream &os); /// binaeres Schreiben
void binRead (Gauss &Gs, ifstream &is); /// binaeres Lesen


ostream& operator << (ostream &os, Gauss gs);

#ifndef OHNE_GAUSS
#define OHNE_GAUSS
const Gauss ohne_gauss={0,Vector<double>(0,0,0),Vector<double>(0,0,0),false};
#endif

/**
    Klasse Strahl:
	repraesentiert einen Strahl mit einem Mittelstrahl und vier Randstrahlen, 
    denen jeweils ein elektrisches Feld E und ein Richtungsvektor k zugeordnet sind.
	Index 4 -> Mittelstrahl
	Ein Strahl kann schrittweise weitergefuehrt werden -> next()
*/

class Ray {

public: 
  Ray();
  Ray(StrahlBuffer &B);
  Ray(Ebene E, const Vector<double> &p, double dy, double dz, const
  Vector<complex<double> > &Pol, complex<double>  n0, double r0, double
  k0, const int Anzein=0, Form **Einschluss=NULL,bool logRay=false);
  Ray(const Vector<double> &p, double dy, double dz, const Vector<complex<double> > &Pol, const Vector<double> &K, complex<double>  n0, double r0,double k0,const int Anzein=0,Form **Einschluss=NULL, bool logRay=false);
  void setGauss (Gauss g) {this->g=g; } 
  bool next(); /// naechster Schritt
  Ray reflect(Vector<double> n, complex<double>  n1, complex<double>  n2); /// Strahl wird reflektiert
  void refract(Vector<double> n, complex<double>  n1, complex<double>  n2); /// Strahl wird gebrochen
  void tunnel(Vector<complex<double> > Pol, complex<double>  n1, complex<double>  n2); 
  void tunnel(Vector<complex<double> > Pol, complex<double>  na, complex<double>  np, int l);
  double flaeche(); /// Flaeche zwischen den Randstrahlen
  void setRefract (complex<double>  n) { this->n=n; } /// Momentanen Brechungsindex setzen
  complex<double>  getRefract() { return n; } /// Gibt Brechungsindex zurueck
  void setk(int i, const Vector<double> &K) { k[i]=K; }
  void setphi(int i,const complex<double>  &p) { phi[i]=p; }
  void setP(int i,const Vector<double> &p) { P[i]=p; }
  void setiR (int i) { iR=i; }
  void setGetunnelt (bool v) { getunnelt=v; }
  void setN0(complex<double> n) {n0=n;}
  Vector<double> getk(int i) { return k[i]; }
  Vector<double> getP(int i) { return P[i]; }
  complex<double>  getphi(int i) { return phi[i]; }
  void getP(Vector<double> p[5]) {p=P;}
  Vector<complex<double> > getAmp(int i)  { return E[i]; }
  void setAmp (int i, const Vector<complex<double> > &A) { E[i]=A; }
  Form *Einschluss(int i) {return Ein[i]; }
  ~Ray();
  bool Einschluss() {return imEinschluss;} 
//  Vector<double> checkEinschluss(int s, const Vector<double>& Ps,int& Index);
void checkEinschluss(int Index[5], Vector<double> *Pmin);

  bool Getunnelt() { return getunnelt; }
  int EinIndex() { return einindex; }
  double Intensity (int i) { return abs(E[i]); }
  friend ostream& operator << (ostream &os,Ray S);
  int reflections () {return iR;}
  double pjump (Vector<double> P1[5],Vector<double> P2[5], const double epsilon=1E-10);
  double pjump (Vector<double> P1[5],Vector<double> P2[5], Vector<double> *S);
  double pjump(void);
  double normVol (Vector<double> P[5],Vector<double> k,  Vector<double> n);
  bool schneideEbene (Vector<double> *Erg, const Ebene &E);
  int schneideEbene(const Ebene &E, double d, Vector<double> &S1, Vector<double> &S2);
  Vector<double> schneideEbene( const Ebene &E, bool &found);
public:
/*!
\param n Oberflaechennormale (in Richtung gegen den Strahl)
\param n1 Brechungsindex Medium 1
\param n2 Brechungsindex Medium 2
Rueckgabewert: gebrochener Strahl
*/
Ray reflect(Vector<double> *n, complex<double>  n1, complex<double>  n2); /// Reflexion des Strahls
/*!
\param n Oberflaechennormale (in Richtung gegen den Strahl)
\param n1 Brechungsindex Medium 1
\param n2 Brechungsindex Medium 2
*/
void refract(Vector<double> *N, complex<double>  n1, complex<double>  n2); /// Brechung des Strahls


  Matrix<complex<double> > Fresnel_reflect (double alpha,complex<double>  n1, complex<double>  n2); 
  Matrix<complex<double> > Fresnel_trans (double alpha, complex<double>  beta, complex<double>  n1, complex<double>  n2);
void init_Efeld (const Vector<complex<double> >& Pol,const int AnzRays, double dx, Ebene Eb); 
void init_EfeldGauss (int i, const Ebene& Eb,
                                   const Vector<complex<double> >& Pol);

  void init_Efeld (const Vector<complex<double> > &Pol,int AnzRays=1);
  Vector<double> crossPlane (const Vector<double> Pe, const Vector<double> n);
  double cross (const Vector<double> P10, const Vector<double> P11,
                       const Vector<double> P20, const Vector<double> P21);
  double crossXAxis (const Vector<double>& P1, const Vector<double>& P2,const Vector<double> &k);
  Vector <double> nextCaustic (double &l);

  Ray& operator = (const Ray& S);
 // project(int i){};
  Vector<double> OK;
  Vector<double> P[5];
  Vector<double> k[5];
  Vector<complex<double> > E[5];
  complex<double>  n,n0;
  double k0;
  double r0,rc;
  complex<double>  phi[5];
  //EinschlussInfo *Ein;
  Form **Ein;
  int AnzEin;
  int einindex;
  int iR;
  int pol;
  bool imEinschluss;
  bool getunnelt;
  bool logRay;
  bool isValid;
  Gauss g;
  Vector<double> ka;
  double KORR;
};


#endif

