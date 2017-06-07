#include "PStrahl.h"


Ray_pow::Ray_pow(void) : IStrahl ()
{
}

Ray_pow::Ray_pow(double pow, const Vector<double> &p,
         const Vector<complex<double> > &Pol, const Vector<double> &K,
         complex<double>  n0, double r0, double k0,
         const int Anzein=0, Form **Einschluss=NULL) : IStrahl(p,Pol, K,n0,r0, k0, Anzein, Einschluss)
{
	this->Pow=pow;
	E1=Pol/abs(Pol);
	E2=E1*sqrt(Pow);
}

Ray_pow Ray_pow::reflect(Vector<double> n, complex<double> n1, complex<double> n2)
/** Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
  * reflektiert (einindex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
{
 Ray_pow Erg;
 Matrix <double> D,H,R;
 Matrix<complex<double> > FR, FT;
 Vector <double> Ph,e0,e1,e2;
 double alpha,gamma;

 // n=-n;

  /* Erst mal die Fresnelmatrix für Reflexion berechnen */

  double nk=n*k/(abs(n)*abs(k));
  double h;
  if (nk<0) { n = -n; nk = -nk; }
  getKSystem(n,k,e0,e1,e2);
  if (nk>1.0) nk=1.0;
  alpha=acos(nk);
  gamma=M_PI-2.0*alpha;
  if (alpha>M_PI/2.0) { gamma=-gamma; alpha=M_PI-alpha;}
  if (einindex>-1) e0=-e0;
  trafo(e0,e1,e2,H,R);
  FR=Fresnel_reflect(alpha,n1,n2);

  /* Fresnelmatrix für Transmission */ 
  FT(0,0)=sqrt(abs(1-abs2(FR(0,0))));
  FT(1,1)=sqrt(abs(1-abs2(FR(1,1))));
  FT(2,2)=sqrt(abs(1-abs2(FR(2,2))));
   D=drehmatrixD(e2,k,gamma);

  
  Erg.E1=this->E1;
  Erg.E2=this->E2;
  Erg.Pow=Pow;

  if (imEinschluss)
  {
  // Strahl will aus dem Einschluss raus
  Erg=*this;
  Erg.OK=zero;
  Erg.imEinschluss=false;
  Erg.n=n2;
  Erg.r0=r0;
  Erg.einindex=-1;  
  Erg.refract(FT,n,n1,n2);
  Erg.Pow=abs2(Erg.E2);
  } // if imEinschluss
  else
  {
   if (einindex!=-1)
   {
    // Einschluss wurde getroffen
    Erg=*this;
    Erg.imEinschluss=true;  
   // Erg.n=Ein[einindex]->n;
    Erg.n=n2;
    Erg.OK=Ein[einindex]->P; // ??
    Erg.einindex=einindex;
	Erg.refract(FT,n,n1,n2);
	Erg.getunnelt=false;
	Erg.Pow=abs2(Erg.E2);
    einindex=-1;
   } // if einindex!=-1
   else // Reflexion an der Partikeloberfläche
   {
     Erg=*this;
     Erg.imEinschluss=false;
     Erg.einindex=-1;
     Erg.n=n2;
     Erg.OK=zero;
     Erg.P=P;
     Erg.E1=E1;
     Erg.E2=E2;
     Erg.k=k;
     Erg.refract(FT,n,n1,n2);
	 Erg.Pow=abs2(Erg.E2);
   }
   /* E1 ist "nur" der Polarisationsvektor und zeigt die Richtung des E-Felds an */

  } // else if imEinschluss
  
   k=D*k;
  iR++;
  E1=D*E1;
  E2=H*E2;
  E2=FR*E2;
  E2=R*E2;
  E2=D*E2;
  
  Pow=abs2(E2);
  
 return Erg; 
}

Ray_pow::~Ray_pow(void)
{
}


void Ray_pow::refract(Matrix<complex<double> > T, Vector<double> N, complex<double> n1, complex<double> n2)
{
 double s;
 Matrix<double> H,R; 
 Matrix<double> D;
 Vector <double> n,e0,e1,e2;
 double alpha, gamma;
 complex<double>  beta;

  n=N/abs(N);
  double nk = (n*k) / abs(k);
  // if (nk > 0) { n = -n;  nk = -nk; }
  getKSystem (n,k,e0,e1,e2);
 
  if (nk>1.0) nk=1.0;
  alpha=acos(nk);
  if (alpha>M_PI/2.0) { alpha=M_PI-alpha; e2=-e2; }
  beta=asin((complex<double>) real(n1)/real(n2)*sin(alpha));
  gamma=real(beta)-alpha;
  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=drehmatrixD(e2,k,gamma);
  
  k=D*k;
  E1=H*E1;
  E1=T*E1;
  E1=R*E1;
  E1=D*E1;
  E1=E1/abs(E1);   // E1 ist nur der Polarisationsvektor und gibt somit nur die Polarisationsrichtung vor, d.h. |E1|=1

  E2=H*E2;
  E2=T*E2;
  E2=R*E2;
  E2=D*E2;
}

ostream & operator << (ostream & os, Ray_pow S)
{
  os << "P=" <<  S.P << "   k=" << S.k << endl;
  os << "E1=" << S.E1 << endl;
  os << "E2=" << S.E2 << endl;
  os << "Pow=" << S.Pow << endl;
  os << "Anzahl Einschlüsse=" << S.AnzEin << endl;
  os << "imEinschluss=" << S.imEinschluss  << "    einindex=" <<
  S.einindex << endl;
  os << "OK=" << S.OK << endl;
 return os;
}