#ifndef SUPERELLIPSOID_N_H
#define SUPERELLIPSOID_N_H

#include "lgsurface.h"
#include "vector.h"
#include "form.h"
#include <math.h>

/**
	@author Thomas Weigel <weigel@lat.rub.de>
*/
class Superellipsoid_n : public LGSurface
{
 /* Vereinfachter Superellipsoid mit 
    (x/a1)^m+(y/a2)^m+(z/a3)^m=0 
*/

public:
    Superellipsoid_n (const Superellipsoid_n &S);
    Superellipsoid_n (const Vector<double> &P, const Vector<double> &r,double m,         
        complex<double>  n,
        double r0,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex=ex,
        const Vector<double> &Ey=ey,
        const Vector<double> &Ez=ez,
        const int type=-1
       );
    Superellipsoid_n();
    ~Superellipsoid_n();
    double getm() {return m;}
    Vector<double> getr() {return r;}    
    void setm (double m) {this->m=m;}
    void setr (const Vector<double> &r) {this->r=r;}
    double f(const Vector<double> &P) ;
    double g(const Vector<double> &P); 
    Vector<double> grad(const Vector<double> &P);
    double calcG(double t1, double t2);
    double geta1() {return r[0]; }
    double geta2() {return r[1]; }
    double geta3() {return r[2]; }
    void setParms (const Vector<double>& r, double m) {this->r=r; this->m=m;}
    void seta1(double v) {r[0]=v;}
    void seta2(double v) {r[1]=v;}
    void seta3(double v) {r[2]=v;}
    void setr0(double r0) 
    {
     r=r/this->r0*r0;
     P=P/this->r0*r0; 
     this->r0=r0;
    } 
  
  
    
    void binWrite (ofstream &os);
    void binRead( ifstream & is);
    void scale(double sf);
    void initQuad();
    double B(double x,double y); // Die Beta-Funktion
    double Volume ();
    

protected :
  double m;
  Vector<double> r;
};

#endif
