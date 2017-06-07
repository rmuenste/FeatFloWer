#ifndef LINSE_H
#define LINSE_H

#define sgn(x) ((x<0) ? -1 : 1)

#include "form.h"
#include <iostream>

class Linse : public Form {
public :
  Linse();
  Linse(complex<double> n, Vector<double> P, Vector<double> N, double r1, double r2, double l1, double l2, double h);
void binWrite (ofstream &os) {}   // MUSS dringend gemacht werden 
void binRead (ifstream &os) {}// MUSS dringend gemacht werden

  void scale(double sf) ; 
  bool next(const Vector<double> &p,const  Vector<double> &k, 
                                 Vector <double> &pout, const  int inside);
Vector<double> norm(const Vector<double> &p);
bool isInside(const Vector<double> &p); 
double Volume () {return -1; }  // muss noch überarbeitet werden !!!!!!!!!
void setr0(double r0);
friend ostream & operator << (ostream &os, Linse); 
void  initQuad() {};
void setP (Vector<double> r); // Setze Bezugspunkt des Einschlusses
  void setP (double x, double y, double z) {setP(Vector<double>(x,y,z)); }
  protected :
  Vector<double> P[2],Pl,N;
  double r[2],l[2],h;
  double r2[2];
};



#endif
