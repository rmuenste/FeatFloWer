/***************************************************************************
 *   Copyright (C) 2005 by Thomas Weigel                                   *
 *   weigel@lat.ruhr-uni-bochum.de                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef ZYLINDER_H
#define ZYLINDER_H

#include "form.h"

/**
Kreiszylinder

	@author Thomas Weigel <weigel@lat.ruhr-uni-bochum.de>
*/
class Zylinder : public Form
{
public:
    Zylinder();
    Zylinder (const Form &F);
    Zylinder (const Zylinder &E);
    Zylinder (
             const Vector<double> &P,
             double r,
             double h,
             complex<double>  n,
             double r0=1.0,
             const Matrix<complex<double> > alpha=CUNITY,
             const Vector<double> &Ex=ex,
             const Vector<double> &Ey=ey,
             const Vector<double> &Ez=ez
             );
    ~Zylinder() {}
    void binWrite (ofstream &os);
    void binRead (ifstream &is);
    void scale (double sf);
    bool next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside=-1);
    Vector<double> norm(const Vector<double> &P);
    bool isInside(const Vector<double> &Ps);
    double Volume();
    void setr0(double r0);
    void setParms(double r, double h);
    void setP (Vector<double> r) {P=r; initQuad();} //
    void setP(double x, double y, double z) {setP(Vector<double> (x,y,z));}
    double radius() {return r/r0;}
    double height() {return h/r0;}
    void setRadius (double r) {this->r=r*r0;}
    void setHeight (double w) {this->h=w*r0;}
    void initQuad();
   double h; // Höhe
   double r; // Radius 
protected:
    double schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP);
    double schneideMantel (Vector<double> &p, Vector<double> &k);
    double schneideDeckel(Vector<double> &p, Vector<double> &k);
};


/** 
  Klasse, die einen Zylinder mit hexagonalem Querschnitt repräsentiert
*/
class ZylinderHexagonal : public Zylinder 
{
 public:
 double a;   
 ZylinderHexagonal(){a=0;h=0;}
 ZylinderHexagonal(double a, double h);
     void binWrite (ofstream &os) {}
    void binRead (ifstream &is) {}
    Vector<double> norm(const Vector<double> &P) { return Vector<double>(); };
    bool isInside(const Vector<double> &Ps) { return false; };
    void scale (double sf){}
    bool next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside = -1) { return false; };
    double Volume() { return 0.0; };
    void setr0(double r0){}
    void setParms(double r, double h){}
    void initQuad(){}
    void setP (Vector<double> r) {P=r; initQuad();} //
    void setP(double x, double y, double z) {setP(Vector<double> (x,y,z));}
friend ostream& operator << (ostream& os, ZylinderHexagonal Z);
// protected:
    Vector<double> e1[6],Pe[6];
    Vector<double> ne[6];
    void init();
//    double schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP);
    double schneideMantel (Vector<double> &p, Vector<double> &k);
//    double schneideDeckel(Vector<double> &p, Vector<double> &k);  
};

#endif
