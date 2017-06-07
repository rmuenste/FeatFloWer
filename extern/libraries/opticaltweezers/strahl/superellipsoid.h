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
#ifndef SUPERELLIPSOID_H
#define SUPERELLIPSOID_H
#include "form.h"
#include "surface.h"


/**
@author Thomas Weigel
*/
class Superellipsoid_D : public surface
{
public:
    Superellipsoid_D();
    Superellipsoid_D (const Superellipsoid_D &Su);
    ~Superellipsoid_D();
Vector<double> super (double theta, double phi,double a1, double a2, double a3, double e1, double e2);
void generate (int ntheta, int nphi, double a1, double a2, double a3, double e1, double e2);
void seta1(double a1); 
void seta2(double a2); 
void seta3(double a3);
void sete1(double e1); 
void sete2(double e2); 
void setntheta(double ntheta); 
void setnphi(double nphi);
void setr0(double r0) {
          this->r0=r0; 
          a1=a1/this->r0*r0; 
	  a2=a2/this->r0*r0; 
	  a3=a3/this->r0*r0; 
            generate(ntheta,nphi,a1,a2,a3,e1,e2); }

double geta1(){return a1;} 
double geta2(){return a2;} 
double geta3(){return a3;}
double gete1(){return e1;} 
double gete2(){return e2;} 
double getntheta(){return ntheta;} 
double getnphi(){return nphi;}
void initQuad() {};
void setP (Vector<double> r) {P=r;initQuad();} // Setze Bezugspunkt des Einschlusses
  void setP (double x, double y, double z) {setP(Vector<double>(x,y,z)); }
void binWrite (ofstream &os);
void binRead (ifstream &is);
// protected :
double a1,a2,a3;
double e1,e2;
double ntheta,nphi;
};


class Superellipsoid : public Form
{
public:
    Superellipsoid();
    Superellipsoid (const Superellipsoid &Su);
Superellipsoid (
             const Vector<double> &P,
             double a1, double a2, double a3, double e1, double e2,
             complex<double>  n,
             double r0=1.0,
             const Matrix<complex<double> > alpha=CUNITY,
             const Vector<double> &Ex=ex,
             const Vector<double> &Ey=ey,
             const Vector<double> &Ez=ez
             );
    ~Superellipsoid(){}
void initQuad();
void binWrite (ofstream &os) {}   // MUSS dringend gemacht werden 
void binRead (ifstream &os) {}// MUSS dringend gemacht werden

double berechneNstNV (double xn, Vector<double> p, Vector<double> k, double epsilon=10e-10);
double berechneG(double &t1, double &t2, Vector<double> p, Vector<double> k);
double berechneg(double &t, Vector<double> p, Vector<double> k);
double berechneFvont(double t, Vector<double> p, Vector<double> k);
double berechneFsvont(double t, Vector<double> p, Vector<double> k);
bool naeherung (double &t1, double &t2, Vector<double> p, Vector<double> k, int &RecCounter);
bool next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside);
double super(double x, double y, double z);
double super(Vector<double> r);
bool isInside (const Vector<double> &p);
double Volume();
Vector<double> norm (const Vector<double> &P); 
void setParms(double a1,double a2, double a3, double e1, double e2);
void scale( double sf);
void setP (Vector<double> r) {P=r;initQuad();} // Setze Bezugspunkt des Einschlusses
  void setP (double x, double y, double z) {setP(Vector<double>(x,y,z)); }
void setr0 (double r0) 
  {
	  a1=a1/this->r0*r0; 
	  a2=a2/this->r0*r0; 
	  a3=a3/this->r0*r0; 
	  this->r0=r0;
  }


double geta1(){return a1;} 
double geta2(){return a2;} 
double geta3(){return a3;}
double gete1(){return e1;} 
double gete2(){return e2;} 

double a1,a2,a3;
double e1,e2;
};
#endif
