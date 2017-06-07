#include "erythrocyte.h"

Erythrocyte::Erythrocyte()
 : LGSurface()
{
 type=ERYTHROCYTE;
 n=1.06;
 ninel=1.06;
r0=1.0;
P=Vector<double>(0,0,0);
 init();
}

void Erythrocyte::init()
{
 d=6.3; // Durchmesser in µm
 a=(0.86*d/2.0)*(0.86*d/2.0);
b=0.01384083;
c=0.2842917;
h=0.01306932;
}
Erythrocyte::~Erythrocyte()
{
}

Erythrocyte::Erythrocyte(const Erythrocyte &E)
:LGSurface (E)
{
 P=E.P;
 n=E.n;
 r0=E.r0;
type=ERYTHROCYTE;
 init();
}

double Erythrocyte::f(const Vector<double> &P) 
{
 double alpha=4.0/(d*d)*(P[0]*P[0]+P[1]*P[1]);
 return P[2]*P[2]-a*(1-alpha)*(b+c*alpha+h*alpha*alpha);
}

 
Vector<double> Erythrocyte::grad(const Vector<double> &P)
{
 double alpha=4.0/(d*d)*(P[0]*P[0]+P[1]*P[1]);
 double A=(2.0*alpha*(a*c-a*h)+3.0*alpha*alpha*a*h+a*b-a*c)*8.0/(d*d);
 return Vector<double> (A*P[0],A*P[1],2.0*P[2]);   
}

double Erythrocyte::calcG(double t1, double t2)
{
 Vector<double> P(Ps+t1*ks);
 double Kx=8.0*ks[0]/(d*d);
 double Ky=8.0*ks[1]/(d*d);
 double alpha=4.0/(d*d)*(P[0]*P[0]+P[1]*P[1]);
 double A=(2*alpha*(a*c-a*h)+3*alpha*alpha*a*h+a*b-a*c)*8.0/(d*d);
 double h1=2.0*(a*c-a*h);
 double h2=64/(d*d*d*d);
 double G1=fabs((Kx+Ky)*A+(h2*(P[0]*P[0]*ks[0]*ks[0]+P[1]*P[1]*ks[1]*ks[1])+Kx*P[0]+Ky*P[1])
                          *(h1+6*a*h*alpha))+2*ks[2]*ks[2]; 
  
 P=Ps+t2*ks;
 alpha=4.0/(d*d)*(P[0]*P[0]+P[1]*P[1]);
 A=(2*alpha*(a*c-a*h)+3*alpha*alpha*a*h+a*b-a*c)*8.0/(d*d);
 double G2=fabs((Kx+Ky)*A+(h2*(P[0]*P[0]*ks[0]*ks[0]+P[1]*P[1]*ks[1]*ks[1])+Kx*P[0]+Ky*P[1])
                          *(h1+6*a*h*alpha))+2*ks[2]*ks[2]; 
 
 if (G1>G2) return G1; 
 return G2;
}

void Erythrocyte::binWrite (ofstream &os)
{
 P.binWrite(os);
 H.binWrite(os);
 R.binWrite(os);
 os.write ((char *) &n, (char) sizeof (n)); 
 alpha.binWrite(os);
 pul.binWrite (os);
 por.binWrite (os);
 for (int i=0; i<3; i++)
  e[i].binWrite(os);
 os.write ((char *) &Ealpha,(char) sizeof (Ealpha));
 os.write ((char *) &Ebeta,(char) sizeof (Ebeta));
 os.write ((char *) &Egamma,(char) sizeof (Egamma));
 os.write ((char *) &sf,(char) sizeof (sf));
 os.write ((char *) &r0,(char) sizeof (r0));
}

void Erythrocyte::binRead (ifstream &is)
{
 type=ERYTHROCYTE;
 P.binRead(is);
 H.binRead(is);
 R.binRead(is);
 is.read((char *) &n, (char) sizeof(n)); 
 alpha.binRead(is);
 pul.binRead (is);
 por.binRead (is);
 for (int i=0; i<3; i++)
  e[i].binRead(is);
 is.read ((char *) &Ealpha,(char) sizeof (Ealpha));
 is.read ((char *) &Ebeta,(char) sizeof (Ebeta));
 is.read ((char *) &Egamma,(char) sizeof (Egamma));
 is.read ((char *) &sf,(char) sizeof (sf));
 is.read ((char *) &r0,(char) sizeof (r0)); 
}