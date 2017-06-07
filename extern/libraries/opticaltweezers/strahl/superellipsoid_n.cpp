#include "superellipsoid_n.h"
#include "goodies.h"
#include <iostream>
using namespace std;

Superellipsoid_n::~Superellipsoid_n()
{
}


Superellipsoid_n::Superellipsoid_n (const Superellipsoid_n &F) : 
LGSurface (F)
{
 r=F.r;
 m=F.m;
 type=SUPERELLIPSOID_N;
}

Superellipsoid_n::Superellipsoid_n() : LGSurface()
{
 type=SUPERELLIPSOID_N;
 r=Vector<double> (0,0,0);
 m=0;
}

 Superellipsoid_n::Superellipsoid_n (const Vector<double> &P, const Vector<double> &r,double m,  
        complex<double>  n,
        double r0,
        Matrix <complex<double> > alpha,
        const Vector<double> &Ex,
        const Vector<double> &Ey,
        const Vector<double> &Ez,
        const int type
       ) : LGSurface(P,n,r0,alpha,Ex,Ey,Ez,SUPERELLIPSOID_N)
{
 this->r=r;
 this->m=m;
 this->r0=r0;
}

double Superellipsoid_n::f(const Vector<double> &P) 
{
 //cout << "pow:" << P << "   r=" << r << "  m" << m << endl;
 return pow(fabs(P[0]/r[0]),m)+pow(fabs(P[1]/r[1]),m)+pow(fabs(P[2]/r[2]),m)-1.0;
 // return P[0]*P[0]/r[0]/r[0]+P[1]*P[1]/r[1]/r[1]+P[2]*P[2]/r[2]/r[2]-1.0;
//  cout << "******" << endl;
}


Vector<double> Superellipsoid_n::grad(const Vector<double> &P)
{
//  cout << "A pow:" << P << "   r=" << r << "  m" << m << endl;
 return /*m*Vector<double> (pow(fabs(P[0]/r[0]),m-1.0)/r[0]*sign(P[0]),
                          pow(fabs(P[1]/r[1]),m-1.0)/r[1]*sign(P[1]),
                          pow(fabs(P[2]/r[2]),m-1.0)/r[2]*sign(P[2]));*/
        2.0*Vector<double> (P[0]/(r[0]*r[0]),P[1]/(r[1]*r[1]),P[2]/(r[2]*r[2]));
// cout << "AAAAA" << endl; 
}


double Superellipsoid_n::calcG(double t1, double t2) 
{
 Vector<double> P1=Ps+ks*t1;
 Vector<double> P2=Ps+ks*t2;
 double G,G1,G2;
 G=0;
 for (int i=0; i<3; i++)
 {
  G1=fabs(m*(m-1)*(ks[i]*ks[i]*pow(fabs(P1[i]/r[i]),m-2.0)/(r[i]*r[i])));
  G2=fabs(m*(m-1)*(ks[i]*ks[i]*pow(fabs(P2[i]/r[i]),m-2.0)/(r[i]*r[i])));
  if (G1>G2) G+=G1; else G+=G2;
 }
//  G=2.0/(r[0]*r[0]);
 return G;
}

void Superellipsoid_n::binWrite (ofstream &os) 
    {      
      P.binWrite(os);
      H.binWrite(os);
      R.binWrite(os);
      os.write ((char *) &n, (char) sizeof (n)); 
      alpha.binWrite(os);
      pul.binWrite (os);
      por.binWrite (os);      
      os.write ((char *) &Ealpha,(char) sizeof (Ealpha));
      os.write ((char *) &Ebeta,(char) sizeof (Ebeta));
      os.write ((char *) &Egamma,(char) sizeof (Egamma));
      os.write ((char *) &sf,(char) sizeof (sf));
      os.write ((char *) &r0,(char) sizeof (r0));
      r.binWrite(os);      
      os.write ((char *) &m, (char) sizeof(m));      
    } 

void Superellipsoid_n::binRead (ifstream &is) 
{
 type=SUPERELLIPSOID_N;
 P.binRead(is);
 H.binRead(is);
 R.binRead(is);
 is.read((char *) &n, (char) sizeof(n)); 
 alpha.binRead(is);
 pul.binRead (is);
 por.binRead (is);
 is.read ((char *) &Ealpha,(char) sizeof (Ealpha));
 is.read ((char *) &Ebeta,(char) sizeof (Ebeta));
 is.read ((char *) &Egamma,(char) sizeof (Egamma));
 is.read ((char *) &sf,(char) sizeof (sf));
 is.read ((char *) &r0,(char) sizeof (r0)); 
 r.binRead(is);      
 is.read ((char *) &m, (char) sizeof(m)); 
}

void Superellipsoid_n::scale(double sf)
{
 this->sf=sf;
 r=r/abs(r)*sf; 
 initQuad();
}

void Superellipsoid_n::initQuad ()
{
   Vector<double> h[8];

   for (int k=0; k<=1; k++)
   for (int l=0; l<=1; l++)
     for (int m=0; m<=1; m++)
     {
        k==0 ? h[k*4+l*2+m][0]=-r[0] : h[k*4+l*2+m][0]=r[0];
        l==0 ? h[k*4+l*2+m][1]=-r[1] : h[k*4+l*2+m][1]=r[1];
        m==0 ? h[k*4+l*2+m][2]=-r[2] : h[k*4+l*2+m][2]=r[2];
     }

 pul=zero;
 por=zero;
 for (int i=0; i<8; i++)
 {
  h[i]=R*h[i];
  if (h[i][0]<pul[0]) pul[0]=h[i][0];
  if (h[i][1]<pul[1]) pul[1]=h[i][1];
  if (h[i][2]<pul[2]) pul[2]=h[i][2];
  if (h[i][0]>por[0]) por[0]=h[i][0];
  if (h[i][1]>por[1]) por[1]=h[i][1];
  if (h[i][2]>por[2]) por[2]=h[i][2];
 }
 pul=pul+P;
 por=por+P;
}

double Superellipsoid_n::B(double x, double y)
{
// return tgamma(x)*tgamma(y)/tgamma(x+y);
	return 0.0; // tgamma muss noch implementiert werden
}

double Superellipsoid_n::Volume ()
{
 double V=2.0*r[0]*r[1]*r[2]*2.0/m*2.0/m*B(1.0/m+1,2.0/m)*B(1.0/m,1.0/m);
 return V;
}