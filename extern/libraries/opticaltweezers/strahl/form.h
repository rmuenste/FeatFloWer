/***************************************************************************
                          form.h  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by Thomas Weigel
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


 
#ifndef FORM_H
#define FORM_H

#ifndef INF
#define INF 1.0/0.0
#endif



#include "vector.h"
#include "matrix.h"
#include "resutil.h"
#include <fstream>
#include <time.h>

/*! Definitionen der unterschiedlichen Objektarten
*/

#define KEINE_FORM  -1
#define ELLIPSOID    0
#define SURFACE      1
#define FUNSURF      2
#define SUPERELLIPSOID_D 3
#define SUPERELLIPSOID   4
#define ZYLINDER     5
#define KREISKEGEL   6
#define KEGELSTUMPF  7
#define COMPOUND     8
#define SPIEGEL      12
#define SUPERELLIPSOID_N 10
#define ERYTHROCYTE      11
#define KEGELSTUMPF_HOHL 9 
#define HOHLFASER    13
#define NINCTYPES    14
#define LINSE	     15
#define ZYLINDER_HEXAGONAL 16
#define EPS 1E-10*r0

using namespace std;
/**
Basisklasse für alle Objekte, die von der Raytracing Bibliothek verwendet wird.
Alle Objekte sind von dieser Klasse abgeleitet.
*@author Thomas Weigel
*/
class Form {
public:
	Form();
  Form (const Form &F);
 // ~Form(){cout << "FORM->DESTRUKTOR" << endl;};  
  /*!
    Hauptkonstruktorvorlage für Form-Objekte:
	P: Ort des Objektes
	n: Brechungsindex 
	alpha: Polarisierbarkeitsmatrix
	Ex,Ey,Ez : Koordinatenachsen des lokalen Koordinatensystems
	type: Art des Objektes
  */
  Form (const Vector<double> &P,                                 
        complex<double>  n,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex=ex,
        const Vector<double> &Ey=ey,
        const Vector<double> &Ez=ez,

        const int type=-1
       );

//  virtual Form& operator = (Form &) = 0;
//  virtual Form *copy () = 0;
  virtual void binWrite (ofstream &os) = 0;                         /// binäres Schreíben in den Filestream os  
  virtual void binRead (ifstream &os) = 0;                          /// binäres Lesen aus dem Filestream os 
  virtual void scale(double sf) = 0;                                /// Skaliert auf ursprüngliche Größe * sf
  virtual bool next(const Vector<double> &p, const Vector<double> &k,
                     Vector<double> &pout,const int inside=-1) = 0; /// Sucht den nächsten Schnittpunkt eines Strahls mit dem Einschluss
  virtual Vector<double> norm (const Vector<double> &P) = 0;        /// Oberflächennormale an der Stelle P
  virtual bool isInside (const Vector<double> &p) = 0;              /// Prüft ob P innerhalb des Einschlusses ist.
  virtual double Volume() = 0;                                      /// gibt Volumen zurück 
  int Type() {return type;}                                         /// gibt Art des Objekttyps zurück
  virtual void initQuad ()= 0;                                      /// Initialisiert den umschriebenen Quader
  virtual  void setr0(double r0)=0;                                 /// Setzt Partikelgröße (damit alle Grössen in Einheiten des Partikelradius angegeben werden können)
  void setMatrix (Matrix<double> H);                                /// Setzt die Matrix zur Transformation ins lokale Koordinatensystem
  void setMatrix (double alpha, double beta, double gamma);         /// Setze Matrizen für Transformation Labor- <-> Einschlusskoordinatensystem 
  virtual void setP (Vector<double> r) = 0; /// Setze Bezugspunkt des Einschlusses
  virtual void setP (double x, double y, double z)= 0; 
  void setn (complex<double> n) {this->n=n; }                       /// setze Brechungsindex
  void setninel (complex<double> ninel)  {this->ninel=ninel; }      /// setze Brechungsindex bei der inelastischen Wellenlänge
  complex<double> getninel () {return ninel; }                      /// gib Brechungsindex bei der inelastischen Wellenlänge zurück   
  complex<double> getn () {return n; }                              /// gib Brechungsindex zurück
  void setPolMatrix (Matrix<complex<double> >alpha) {this->alpha=alpha;}   /// setze Polarisierbarkeitsmatrix
  bool isActive () { return inelactive;}                           
  void setActive (bool active) {inelactive=active;}
  void setAlpha(double Alpha) {setMatrix( Alpha,Ebeta,Egamma); }   /// Setze Drehwinkel um x-Achse (Alpha)
  void setBeta(double Beta) {setMatrix( Ealpha,Beta,Egamma); }     /// Setze Drehwinkel um y-Achse (Beta) 
  void setGamma(double Gamma) {setMatrix( Ealpha,Ebeta,Gamma); }   /// Setze Drehwinkel um z-Achse (Gamma)
  Vector<double> P;                  ///< Ort (=Ursprung des Koordinatensystems)
  Matrix<double> H,R;                ///< Matrizen zur Umrechnung  Einschluss <-> Laborsystem
  complex<double>  n;                ///< Brechungindex (einfallendes Feld);
  complex<double> ninel;             ///< Brechungindex (inelastisch);
  Matrix<complex<double> > alpha;    ///< Polarisierbarkeitsmatrix
  int type;                          ///< Typ des Einschlusses
  Vector<double> pul,por; /// Ecke unten(oben) vorne(hinten) links(rechts) 
  Vector<double> e[3];    /// Einheitsvektoren
  double Ealpha,Ebeta,Egamma; /// Winkel, um die das Objekt gedreht wurde (um x-, dann um die y- und schliesslich um die z-Achse
  double r0;                  /// Radius der "Weltkugel", die das Berechnungsgebiet vorgibt
  double sf;         /// Skalierungsfaktor     
  bool inelactive;   /// Ist das Objekt inelastisch aktiv ?
};

Matrix<double> computeInertia(Form *F); /// Berechne Trägheitsmatrix
// #include "misc.h"

#endif
