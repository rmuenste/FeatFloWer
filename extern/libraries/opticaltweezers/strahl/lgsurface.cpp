#include "lgsurface.h"



LGSurface::LGSurface (const Vector<double> &P,
        complex<double>  n,
        double r0,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex,
        const Vector<double> &Ey,
        const Vector<double> &Ez,
        const int type) :Form (P,n,alpha,Ex,Ey,Ez,type)
{
}

LGSurface::LGSurface (const LGSurface &F) : Form ( F)
{
 type=-1;
}

LGSurface::LGSurface() : Form()
{
 type=-1;
}

double LGSurface::g(const Vector<double> &p)
{
 return ks*grad(p-P);
}

 

LGSurface::~LGSurface()
{
}

bool LGSurface::hatNS (double &a, double &b, int &c)
{
 if ((c>100000)|| ((b-a)/b<1E-10)) 
 {  
   if (f(a)*f(b)<0) return true;
   else return false;
  return true;
 }
 c++;
 
 double ah,bh;
 double tm=(a+b)/2.0;
 double d=(b-a)/2.0;
 double G;
//  cout << "ps=" << Ps << "   ks=" << ks << "   f(Ps)=" << f(Ps) << endl;
// cout << "hatNS: c=" << c << "   a=" << a << "   b=" << b << " f(a)=" << f(a) << "   f(b)=" << f(b) << endl;
 
 double fa=f(a);
 
 double fb=f(b);
 

 if (fa*fb<0)
 { 
  G=calcG(a,b);
  if (fabs(g(tm))<G*d)  // g hat keine NS
     {  
       return true; 
      }
  else
   {
    ah=a; bh=b; b=tm;
    if (hatNS(a,b,c)) return true;
    else 
    {
      a=tm; b=bh;
      if (hatNS(a,b,c)) return true;
      else return false;
     }
   }  
 }
 else
 { 
   G=calcG(a,b);
 if (fabs(g(tm))>G*d)  // g hat keine NS
     {  
       return false; 
      }
  else
   {    
    ah=a; bh=b; b=tm;
    if (hatNS(a,b,c)) return true;
    else 
    {
      a=tm; b=bh; 
      if (hatNS(a,b,c)) return true;
      else return false;
     }
   }   
  return false;
 } 
}


bool LGSurface::next(const Vector<double> &p, const Vector<double> &k,
                     Vector<double> &pout,const int insidE)
{
 double t;
 double a,b=2.0*r0;
 
 ks=H*k;
 Ps=H*(p-P);
 a=0;
 if (fabs(f(a))<1E-10)
 {
  a=1E-8/fabs(g(a))+a;
 }

 int c=0;
 if (hatNS(a,b,c)) 
 {
  t=berechneNS(a,b);
  pout=p+t*k;
  return true;
 }
 else
 {
  pout=p;
  return false;
 }
}


double LGSurface::berechneNS(double a, double b)
{
 
 // Bisection 
 double c,fa;
 do
 {  
  c=(a+b)/2.0;
  if (f(c)*f(a)>0) a=c;
  else b=c; 
 } while (fabs((a-b)/b)>1E-10);
 return c;
 
 /*// Newton
 double x0,x1;
 x0=(a+b)/2.0;
 do
  x1=x0-f(x0)/g(x0);
 while (fabs(x0-x1)>1E-10);
 return x1;*/ 
}
