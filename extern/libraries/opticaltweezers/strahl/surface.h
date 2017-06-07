#ifndef SURFACE_H
#define SURFACE_H

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "vector.h"
#include "form.h"
#include <string.h>

#ifndef SQRT3
#define SQRT3 1.73205080756887729352744634150587236694280525381038062805580  /// SQRT(3)
#endif

/** \page refdreieck
  Klasse, die ein Dreieck repräsentiert. Wird in der Klasse \ref refsurface "surface". */
class dreieck
{


public:  
//  protected:
	  double u, v;
  Vector<double> P[3];  ///< Eckpunkte des Dreiecks
  Vector<double> n;     ///< Oberflächennormale
  public:
	  double area()
	  {
		  Vector<double> h1, h2;
		  h1 = P[1] - P[0];
		  h2 = P[2] - P[0];
		  return abs(h1%h2) / 2.0;
	  }
  int schnittpunkt(Vector<double> r, Vector<double> k, Vector<double> &p, double  eps=1.e-8);
  int schnittpunkt(Vector<double> r, Vector<double> k, double &t, Vector<double> &p, double  eps=1.e-8);
  void setnorm(void);   ///< Berechnet die Normale aus den Eckpunklten  
  void setnorm(Vector<double> n) {this->n=n;} ///< Setzt die Oberflächennormale auf n
  Vector<double>  getnorm(void); ///< gibt die Oberflächennormale zurück
  dreieck();
  dreieck(Vector<double> P1, Vector<double> P2, Vector<double> P3);
  dreieck(Vector<double> ip1, Vector<double> ip2, Vector<double> ip3, Vector<double> P);
  Vector<double>& operator[](int i);
  const Vector<double>& operator[](int i) const;
  dreieck& operator=(const dreieck &dr);
  void binWrite (ofstream &os);
  void binRead (ifstream &is);
  ~dreieck(); 
  friend class surface;
};

ostream& operator << (ostream &os, const dreieck &dr);


/** Klasse, die eine Form repräsentiert, die aus Dreiecken zusammengebaut ist. Die Dreiecke werden durch die Klasse \ref refdreieck "dreieck" 
beschrieben */
class surface : public Form
{
  public:

// Konstruktoren und Destruktor
  surface(int anz);
  surface();
  surface(Vector<double> Oh);
  surface(const Form &);
  surface(const surface &Su);
  surface(const Vector<double> &P, 
                complex<double>  n,
          const Matrix<complex<double> > alpha=CUNITY,
          const Vector<double> &Ex=ex,
          const Vector<double> &Ey=ey,
          const Vector<double> &Ez=ez
         );
  surface(const Vector<double> &P, 
                complex<double>  n,
                int anz, dreieck* list,
          const Matrix<complex<double> > alpha=CUNITY,
          const Vector<double> &Ex=ex,
          const Vector<double> &Ey=ey,
          const Vector<double> &Ez=ez);

  ~surface();
  Vector<double> calcCoM(); /// Schwerpunkt berechnen
  dreieck* S; /// Liste der Dreiecke
  int anzp; /// Anzahl Dreiecke

//   double a;

// Erzeugen der Dreiecksliste
// 1. interaktiv
  int createsurface();
  /** Skaliere die Dreiecke */
  void scale (double sf);

// 2. Einlesen aus Datei FName
  int createsurface(char* FName); /// lädt eine Datei mit dem Namen "FName"
  int importBinSTL(char *FName); ///< importiere binäre STL Datei
 void exportSRF (char *FName);   ///< exportiere als SRF-Datei
// Andere Hilfsfunktionen
/// 1. Rückgabe einer Klasse mit leerem S
  surface nosurface();
/// 2. Hinzufügen einer Liste S
  void addS(dreieck *S, int anz);
/// 2. Löschen von S
  void clearS();
//  string getFName() {return FName; }
  char *getFName() {return FName; }
  void setFilename(char *FName)
  {
    if (this->FName!=0) delete[] this->FName;
    this->FName=new char[strlen(FName)+1];
    strcpy(this->FName,FName);
  }
// Operatoren 
  surface& operator=(const surface&);

// Die virtuellen Funktionen der Form-Klasse
  bool next(const Vector<double> &r, const Vector<double> &k, Vector<double>& p, const int inside=-1);
  Vector<double> norm (const Vector<double> &P);
  bool isInside (const Vector<double> &p);
  void setr0(double r0);
  void initQuad();
  Matrix<double> computeInertia();
void setP (Vector<double> r); // Setze Bezugspunkt des Einschlusses
  void setP (double x, double y, double z) {setP(Vector<double>(x,y,z)); }
  /** No descriptions */
  double isInHost(void);
  void binWrite (ofstream &os);
  void binRead (ifstream &is);
  double Volume();
  double calcVolume();
  void setCenter(Vector<double> P);
  void setCenter2CoM();
  int getCurrentIndex() { return currentIndex; }
  dreieck& getTriangle(int i) { return S[i]; }

  protected:
  Vector<double> currentnorm;  
  int currentIndex;
  char *FName;
};

ostream& operator << (ostream &os, const surface &su);

dreieck operator + (const dreieck &dr, const Vector<double> &v);

dreieck operator - (const dreieck &dr, const Vector<double> &v);

dreieck operator + (const Vector<double> &v, const dreieck &dr);

dreieck operator - (const Vector<double> &v, const dreieck &dr);

dreieck operator * (const Matrix<double> &M, const dreieck &dr);

dreieck operator * (const dreieck &dr, double a);

dreieck operator * (double a, const dreieck &dr);

dreieck operator / (const dreieck &dr, double a);

surface operator + (const surface &s, const Vector<double> &v);

surface operator - (const surface &s, const Vector<double> &v);

surface operator + (const Vector<double> &v, const surface &s);

surface operator - (const Vector<double> &v, const surface &s);

surface operator * (const Matrix<double> &M, const surface &s);

surface generateHexagonCylinder(double a, double h, Matrix<double> M=UNITY); /// Generiert einen hexagonalen Zylinder mit Seitenlänge a und Höhe h


#endif
