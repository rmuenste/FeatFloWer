#ifndef LGSURFACE_H
#define LGSURFACE_H

#include <form.h>

/**
	@author Thomas Weigel <weigel@lat.rub.de>
*/
class LGSurface : public Form
{
public:
     LGSurface ();
     LGSurface (const LGSurface &F);
     LGSurface (const Vector<double> &P,
        complex<double>  n,
        double r0,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex=ex,
        const Vector<double> &Ey=ey,
        const Vector<double> &Ez=ez,
        const int type=-1
       );
    ~LGSurface();
    virtual double calcG(double t1, double t2) =0;
    virtual double f(const Vector<double> &P) = 0;
    virtual Vector<double> grad(const Vector<double> &P) = 0; 
    double g(const Vector<double> &); 
    Vector<double> norm (const Vector<double> &p) 
    {
      Vector<double> n=R*grad(H*(p-P));
      return n/abs(n);
    } 

    
    
    bool next(const Vector<double> &p, const Vector<double> &k,
                     Vector<double> &pout,const int inside=-1);
    bool isInside (const Vector<double> &P) {return f(P)<=0;}
    bool hatNS (double &a, double &b, int &c);
    void setP(Vector<double> r) {P=r; }
    void setP(double x, double y, double z) {setP(Vector<double>(x,y,z));}
    double berechneNS(double a, double b);
  void setr0(double r0) {P=P/this->r0*r0; this->r0=r0; }
 protected :
  double f(double t) {return f(Ps+t*ks);};  // implizite Funktion
  double g(double t) {return g(Ps+t*ks);};  // g(t)=df/dt
  Vector<double> Ps,ks; // Ort und Richtung des Strahls (zur Zwischenspeicherung)
  void binWrite( ofstream &os) {}
  void binRead( ifstream &ss) {}
  double Volume() {return 0;}
  void scale (double sf) {}
  void initQuad() {}
 
  
};

#endif
