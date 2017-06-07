#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <complex>
#include "vector.h"
#include "matrix.h"
#include "form.h"
#include <iostream>
#include <fstream>
/*!
  Von Form abgeleitete Klasse. Repräsentiert ein Ellipsoid
  Die Oberfläche wird beschrieben durch: \f$\frac{x^2}{a^2}+\frac{y^2}{b^2}+\frac{z^2}{c^2}=1\f$
*/

class FormEllipsoid : public Form
{
 public :
 // FormEllipsoid &operator = (const FormEllipsoid& e);
  ~FormEllipsoid(){}; 
  FormEllipsoid ();
  FormEllipsoid (const Form &);
  FormEllipsoid (const FormEllipsoid &E);
  FormEllipsoid (
             const Vector<double> &P,
             const Vector<double> &r,
             complex<double>  n,
             double r0=1.0,
             const Matrix<complex<double> > alpha=CUNITY,
             const Vector<double> &Ex=ex,
             const Vector<double> &Ey=ey,
             const Vector<double> &Ez=ez
             );
  void scale(double sf); /// Skaliere Objekt
  bool next(const Vector<double> &Ps, const Vector<double> &K,
                     Vector<double> &pout,const int inside=-1); /// Berechne nächsten Schnittpunkt mit dem Strahl, repräsentiert durch Ps (Ort) und k (Richtung)
  FormEllipsoid & operator =  (FormEllipsoid &f);
  FormEllipsoid & operator = (FormEllipsoid f);
//  void setMatrix (double alpha, double beta, double gamma);
  Vector<double> norm (const Vector<double> &p); /// Oberflächennormale am Punkt p
  bool isInside (const Vector<double> &p);
//ir  Vector<double> r;   // "Radien der Ellipse" r=(a,b,c) fuer: (x/a)^2+(y/b)^2+(z/c)^2=1
  FormEllipsoid * copy ();
  double Volume ();
  void copy(FormEllipsoid *E);
  void initQuad ();
  friend ostream &operator << (ostream &os,FormEllipsoid E);
 void binWrite(ofstream &os);
  void binRead(ifstream &os); 
  double a() { return r[0]; }
  double b() { return r[1]; }
  double c() { return r[2]; }
  Vector<double> getr () {return r;}  
  void setP(Vector<double> r) {P=r; }
  void setP(double x, double y, double z) {setP(Vector<double>(x,y,z));}
  void setr (Vector<double> &r);
  void setr (double a,double b, double c);
  void seta(double a, bool VConst=false ); 
  void setb(double b, bool VConst=false );
  void setc(double c, bool VConst=false );
  void setr0(double r0);
   Matrix<double> computeInertia ()   /// Berechne Trägheitsmatrix
   {
        Matrix<double> I;
        I(0,0)=1.0/5.0*(r[1]*r[1]+r[2]*r[2]);
        I(1,1)=1.0/5.0*(r[0]*r[0]+r[2]*r[2]);
        I(2,2)=1.0/5.0*(r[0]*r[0]+r[1]*r[1]);
        return I / 1E-12; // 1E-12, da r hier in µm angegeben wird, ich will aber I in m²
   }
// protected :
 Vector<double> r;
 Vector<double> r_2;      // r_2=(1/r[0]^2,1/r[1]^2,1/r[2]^2)
 Vector<double> P2;      //  P2=(P[0]^2,P[1]^2,P[2]^2)
};
#endif
