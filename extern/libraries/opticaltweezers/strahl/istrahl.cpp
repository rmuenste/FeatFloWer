/***************************************************************************
                          strahl.cpp  -  description
                             -------------------
    begin                : Fri Oct 15 1999
    copyright            : (C) 1999 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "istrahl.h"
extern Matrix<double> Dall;
IStrahl::IStrahl()
{
 iR=0;
 einindex=-1;
 OK=zero;
}

IStrahl::~IStrahl(){
}

void IStrahl::checkEinschluss(int &Index, Vector<double> &Pmin)
{
 Vector<double> P1,n;
 double amin,a;
 bool found;
 bool isInside=(einindex>-1);
  amin=100.0*r0;
  Pmin=P;
  Index=-1;
  for (int i=0; i<AnzEin; i++)
  {
   found=Ein[i]->next(P,k,P1,isInside);
   if (found)
   {
    a=abs(P1-P);
    if (a<amin)
    {
      amin=a;
      Pmin=P1;
      /*if (i!=einindex)*/ Index=i;
    }
   }
  }
  
 
   /*found=true;
   for (int j=0; (j<=3) && found; j++)
    for (int i=j+1; (i<=4) && found; i++)
     found=found && (Index[j]==Index[i]);*/
//    isValid=found; 

 } 




/*Vector<double> IStrahl::checkEinschluss(const Vector<double>& P,int& index)
{
 Vector<double> P1,Pmin;
 double amin,a;
 bool found;

 amin=1E+99;
 index=-1;
 Pmin=P;
 for (int i=0; i<AnzEin; i++)
 {
 //  P1=nextP(P, k, Ein[i].P*r0,Ein[i].a*r0,found);
   found=Ein[i]->next(P,k,P1);
   if (found)
   {
    a=abs(P1-P);
    if (a<amin)
    {
      amin=a;
      Pmin=P1;
      if (i!=einindex) index=i;
    }
   }
  }
 return Pmin;
}*/

IStrahl::IStrahl(const Vector<double> &p,
         const Vector<complex<double> > &Pol, const Vector<double> &K,
         complex<double>  n0, double r0, double k0,
         const int Anzein, Form **Einschluss)
{
 double l;
 E1=czero;
 E2=czero;
 P=p;
 k=K;
 Ein=Einschluss;
 AnzEin=Anzein;
 n=n0;
 this->r0=r0;
 this->k0=k0;
 l=k0/(2.0*M_PI);
// KORR=1E-10;
 //init_Efeld(Pol);
 imEinschluss=false;
 einindex=-1;
 OK=zero;
 getunnelt=false;
}


void IStrahl::next()
{ 
 int Index;
 Vector<double> p,R;  
 bool found=true;
 bool Einschluss=true;
 if (!imEinschluss) // Strahl ist ausserhalb eines Einschlusses
 {
   if (AnzEin==0) Index=-1;  
   else checkEinschluss(Index,R);  // Schnittpunkt mit Einschluss suchen
  if(Index==-1)    // Es wurde kein Einschluss gefunden !
  {   
    R=nextP(this->P,k,zero,r0,found); // Schnittpunkt mit der Außenkugel
	if (!found) { R=this->P; Index=-2; } // Wenn nicht gefunden => Da stimmt was nicht 
  } 
  
  einindex=Index;  
    E1=E1*exp(I*k0*n*abs(R-P));   
    E2=E2*exp(I*k0*n*abs(R-P));
  //  this->P=R;
 }
 else
 {   
   Ein[einindex]->next(P,k,R,true);
    E1=E1*exp(I*k0*Ein[einindex]->n*abs(R-P));
    E2=E2*exp(I*k0*Ein[einindex]->n*abs(R-P));
 }

 this->P=R;
 /*if (abs(this->P-R)>1E-15) this->P=R;
 else einindex=-3;*/
}

/*
void IStrahl::next()

//     Hier wird der Strahl zum naechsten Punkt in Strahlrichtung gebracht.
//     Dabei wird ueberprueft,ob ein Einschluss getroffen wurde

{
  int Index;
  Vector<double> R,p,PE;
  bool found=true;
    if (!imEinschluss)
 {
  einindex=-1;
   Index=-1;
   p=nextP(P, k, zero, r0,found);
   R=p;
   PE=checkEinschluss (P,Index);
   if ((Index>-1) && ( abs(PE-P)<abs(p-P)))
     {
      R=PE;
     }
    E1=E1*exp(I*k0*n*abs(R-P));
    E2=E2*exp(I*k0*n*abs(R-P));
    P=R;
 } // if
 else
 {
   Ein[einindex]->next(P,k,R);
   E1=E1*exp(I*k0*Ein[einindex]->n*abs(R-P));
   E2=E2*exp(I*k0*Ein[einindex]->n*abs(R-P));
   P=R;
 }
 if (!imEinschluss) einindex=Index;
}
*/
void IStrahl::refract(Vector<double> N, complex<double>  n1, complex<double>  n2)
{
 double s;
 Matrix<double> H,R;
 Matrix<complex<double> > T;
 Matrix<double> D;
 Vector <double> n,e0,e1,e2;
 double alpha, gamma;
 complex<double>  beta;

  n=N/abs(N);
  getKSystem (n,k,e0,e1,e2);
  double nk=(n*k)/abs(k);
   alpha=acos(nk);
  if (alpha>M_PI/2.0) { alpha=M_PI-alpha; e2=-e2; }
  beta=asin((complex<double>) real(n1)/real(n2)*sin(alpha));
 gamma=real(beta)-alpha;
  //if (real(n2)<real(n1)) { e2=-e2;}
  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=drehmatrixD(e2,k,gamma);
  if (beta!=0.0) T=Fresnel_trans(alpha,beta,n1,n2);
  else
  {
   T(0,0)=(2.0*n1)/(n2+n1);
   T(1,1)=T(0,0);
   T(2,2)=T(0,0);
  }

  k=D*k;
  E1=H*E1;
  E1=T*E1;
  E1=R*E1;
  E1=D*E1;

  E2=H*E2;
  E2=T*E2;
  E2=R*E2;
  E2=D*E2;
}

void IStrahl::tunnel(Vector<complex<double> > Pol, complex<double>  n1, complex<double>  n2)
{
/* double s;
 Matrix<double> H,R;
 Matrix<complex<double> > T;
 Matrix<double> D;




  E=Pol;
  b=sqrt(P[1]*P[1]+P[2]*P[2]);
  n=P/abs(P);
  getKSystem (n,k,e0,e1,e2);
  alpha=acos(n*k);
  if (alpha>M_PI/2.0)
  {
   alpha=M_PI-alpha;
   beta=asin(n1/n2*sin(alpha));
   e2=-e2;
  }

  beta=asin(n1/n2*sin(alpha));
  E*=exp(I*n1*k0*P[0]);
  P[0]=0.0;
  P*=r0/abs(P);
  gamma=-real(acos(b*n1/(r0*n2)));

  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=drehmatrixD(e2,k,gamma);
  T=Fresnel_trans(alpha,beta);
  T*=exp(I*k0*cos(beta)*(b-r0)*n1);
  k=D*k;
  E=H*E;
  E=T*E;
  E=R*E;
  E=D*E;
  P*=(1.0-KORR);*/
}

void IStrahl::tunnel(Vector<complex<double> > Pol, complex<double>  np, complex<double>  na, int l)
{
 double s, gamma;
 Matrix<double> H,R;
 Matrix<complex<double> > T;
 Vector<double> n, e0, e1, e2;
 Matrix<double> D;
 
 // Innerer Reflexionswinkel (von der Tangente an P aus gerechnet)
// cout << "l:" << l << ", np:" << np << ", k0:" << k0 << ", r0:" << r0 << "\n";
 gamma = acos((l+0.5)/(real(np)*real(k0)*r0));
//  cout << "gamma:" << gamma << "\n";
   
 // Rechtsrum oder linksrum tunneln ?
   
 // noch zu implementieren
    
 // P liegt auf der Kugeloberflaeche
 n = P/abs(P);
// cout << "n:" << n << "\n";

 // Transformation in lokales Koordinatensystem
 getKSystem (n,k,e0,e1,e2);
 trafo(e0,e1,e2,H,R);
 // Drehmatrix D
 D=drehmatrixD(e2,k,gamma);
 // k drehen
 k=D*k; 

  E1=H*E1;
//  E1=T*E1;
  E1=R*E1;
  E1=D*E1;

  E2=H*E2;
//  E2=T*E2;
  E2=R*E2;
  E2=D*E2;

/*  E=Pol;
  b=sqrt(P[1]*P[1]+P[2]*P[2]);
  n=P/abs(P);
  getKSystem (n,k,e0,e1,e2);
  alpha=acos(n*k);
  if (alpha>M_PI/2.0)
  {
   alpha=M_PI-alpha;
   beta=asin(n1/n2*sin(alpha));
   e2=-e2;
  }

  beta=asin(n1/n2*sin(alpha));
  E*=exp(I*n1*k0*P[0]);
  P[0]=0.0;
  P*=r0/abs(P);
  gamma=-real(acos(b*n1/(r0*n2)));

  s=1.0;
  trafo(e0,e1,e2,H,R);
  D=drehmatrixD(e2,k,gamma);
  T=Fresnel_trans(alpha,beta);
  T*=exp(I*k0*cos(beta)*(b-r0)*n1);
  k=D*k;
  E=H*E;
  E=T*E;
  E=R*E;
  E=D*E;
  P*=(1.0-KORR);*/
}

IStrahl IStrahl::reflect(Vector<double> n, complex<double>  n1, complex<double>  n2)
/* Strahl wird an einer Oberflaeche reflektiert. Wird an einer Einschlussoberflaeche
   reflektiert (einindex >-1 ), dann wird der transmittierte Strahl zurueckgegeben */
{
 IStrahl Erg;
 Matrix <double> D,H,R;
 Matrix<complex<double> > F;
 Vector <double> Ph,e0,e1,e2;
 double alpha,gamma;
  
  Erg.E1=this->E1;
  Erg.E2=this->E2;

  if (imEinschluss)
  {
  // Strahl will aus dem Einschluss raus
  Erg=*this;
  Erg.OK=zero;
  Erg.imEinschluss=false;
  Erg.n=n2;
  Erg.r0=r0;
  Erg.einindex=-1;
  Erg.refract(n,n1,n2);
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
    Erg.OK=Ein[einindex]->P;
    Erg.einindex=einindex;
	 Erg.refract(n,n1,n2);
	Erg.getunnelt=false;
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
     Erg.refract(n,n1,n2);
   }
  } // else if imEinschluss

  n=-n;
  double nk=n*k/(abs(n)*abs(k));
  getKSystem(n,k,e0,e1,e2);
 // cout << "nk=" << nk << "   n=" << n << "   k=" << k << "   P=" << P/abs(P) << endl; 
  /*if (fabs(nk)>=1) { alpha=0.0;  }
  else */  alpha=acos(nk);
  gamma=M_PI-2.0*alpha;
  if (alpha>M_PI/2.0) { gamma=-gamma; alpha=M_PI-alpha;}
  if (einindex>-1) e0=-e0;
  trafo(e0,e1,e2,H,R);
  F=Fresnel_reflect(alpha,n1,n2);
  D=drehmatrixD(e2,k,gamma);

  k=D*k;

  E1=H*E1;
  E1=F*E1;
  E1=R*E1;
  E1=D*E1;

  E2=H*E2;
  E2=F*E2;
  E2=R*E2;
  E2=D*E2;

  iR++;
 return Erg; 
}

Matrix<complex<double> > IStrahl::Fresnel_reflect (double alpha, complex<double>  n1, complex<double>  n2)
{
 complex<double>  rs,rp;
 Matrix<complex<double> > R;
 double  n12;
 double  beta;

 if (alpha!=0.0)
 {
 n12=real(n2)/real(n1);
 beta=real(asin((complex<double> ) sin(alpha) / n12));
 rp=tan(alpha-beta)/tan(alpha+beta);
 rs=-sin(alpha-beta)/sin(alpha+beta);
 }
 else
 {
  rp=(real(n2-n1))/real(n1+n2);
  rs=-rp;
 }

 R(0,0)=rp;
  R(1,1)=rp;
  R(2,2)=rs;
 return R;
}

Matrix<complex<double> > IStrahl::Fresnel_trans (double alpha,complex<double>  beta, complex<double>  n1, complex<double>  n2)
{
 complex<double>  ts,tp;
 Matrix<complex<double> > T;
 /*Matrix<double> H,R;
 Vector<double> e1,e2,e3;
 complex<double>  n12;
 double alpha;
 complex<double>  beta,Erg;
 alpha=acos(k*n/(abs(k)*abs(n)));
 if (alpha>M_PI/2.0) alpha=M_PI-alpha;
 n12=n2/n1;
 beta=asin((complex<double> ) sin(alpha) / n12);*/

 if (alpha==0.0)
 {
  ts=2.0*real(n1)/real(n1+n2);
  tp=ts; 
 }
 else 
 {
  ts=2.0 * sin(beta) * cos(alpha) / sin(alpha+beta);
  tp=2.0 * sin(beta) * cos(alpha) / (sin(alpha+beta)*cos(alpha-beta));
 }

 T(0,0)=tp;
 T(1,1)=tp;
 T(2,2)=ts;
 
 return T;
}

void IStrahl::init_Efeld (const Ebene& Eb, const Vector<complex<double> >& Pol,const int AnzRays)
{
  double x1,x2,xn,p;
  Vector<double> r;
  k/=abs(k);
  x1=P*Eb.e1;
  x2=P*Eb.e2;
  xn=P*Eb.n;
  r=P;
  p=x1*x1+x2*x2;
  P=x1*Eb.e1+x2*Eb.e2-sqrt(fabs(r0*r0-p))*Eb.n;  
   E1=Pol*exp(I*abs(r-P)*k0);
   E2=Pol*exp(I*abs(r-P)*k0);
}

void IStrahl::init_Efeld (const Ebene& Eb, const Vector<complex<double> >& Pol1, const Vector<complex<double> >& Pol2, const int AnzRays)
{
  double x1,x2,xn,p;
  Vector<double> r;
  k/=abs(k);
  x1=P*Eb.e1;
  x2=P*Eb.e2;
  xn=P*Eb.n;
  r=P;
  p=x1*x1+x2*x2;
  P=x1*Eb.e1+x2*Eb.e2-sqrt(fabs(r0*r0-p))*Eb.n;  
   E1=Pol1*exp(I*abs(r-P)*k0);
   E2=Pol2*exp(I*abs(r-P)*k0);
}



void IStrahl::init_Efeld (
                                   const Vector<complex<double> >& PolS,
                                   const Vector<complex<double> > &PolP,
                                   const int AnzRays)
{
  Vector<double> hr,h;
  hr=P;
  h=P-(P*k)*k;
  P=h-sqrt(r0*r0-abs(h)*abs(h))*k;
  E1=PolS*exp(I*abs(hr-P)*k0); 
  E2=PolP*exp(I*abs(hr-P)*k0); 
}

void IStrahl::init_EfeldGauss (const Ebene& Eb,
                                   const Vector<complex<double> >& PolS,
                                   const Vector<complex<double> > &PolP,
                                   Gauss g)
{ 
 k=g.F-P;
 k=k/abs(k);
 double theta;
 double z=abs(P-g.F);
 double r;
 double R;
 double l0=2.0*M_PI/real(k0);
  Vector<double> hr,h;
  complex<double> E0;
  double w,G;
 
  double z0=M_PI*g.w0*g.w0/l0;
  theta=atan2(g.w0,z0);
  r=sin(theta/2.0)*z;
  w=g.w0*sqrt(1.0+z*z/(z0*z0));
  R=z*sqrt(1.0+z*z/(z0*z0));
  G=atan2(z,z0);
  E0=g.w0/w*exp((-r*r/w*w)-I*k0*r*r/(2*R)-I*k0*z+I*G);  
  double l;
  double Pk=P*k;
  l=-Pk-sqrt(Pk*Pk-P*P+r0*r0);
  P=P+l*k;  
  E1=PolS*E0*exp(I*k0*l);
  E2=PolP*E0*exp(I*k0*l);
}


 double IStrahl::cross (const Vector<double> P10, const Vector<double> P11,
                       const Vector<double> P20, const Vector<double> P21)
 {
  Vector<double> a0,a1;
  double k1,k2,k;

  a0=P20-P10;
  a1=P21-P11;
  k1=a0[1]/a0[0]*(P11[0]-P10[0])+P10[1]-P11[1];
  k2=a1[1]-a0[1]/a0[0]*a1[0];
  k=k1/k2;
  return k;
 }

 double IStrahl::crossXAxis (const Vector<double>& P1, const Vector<double>& P2,
   const Vector<double>& k)
 {
  double m=0.0;

  if (k[0]!=0.0) m=-P1[0]/k[0];
  else m=-1000000.0;
  return m;
 }

Vector<double> IStrahl::crossPlane (const Vector<double> Pe, const Vector<double> n)
/* Berechnet den Schnittpunkt des zentralen Strahls mit einer Ebene, die durch den Aufpunkt Pe und die Normale n beschrieben wird 
 * Die Routine liefert inf-Vektor, falls kein Schnittpunkt */
{
 double l=k*n;
 double h=(Pe-P)*n;
 Vector<double> S;

 if ( (l==0) && (h==0) ) return P; // Strahl befindet sich in der Ebene und bewegt sich in der Ebene
 if ( (l==0) && (h!=0) ) return Vector<double> (-1,-1,-1);
 l=h/l;
 S=P+l*k;
 return S;
}


ostream & operator << (ostream & os, IStrahl S)
{
  os << "P=" <<  S.P << "   k=" << S.k << endl;
  os << "E1=" << S.E1 << endl;
  os << "E2=" << S.E2 << endl;
  os << "Anzahl Einschlüsse=" << S.AnzEin << endl;
  os << "imEinschluss=" << S.imEinschluss  << "    einindex=" <<
  S.einindex << endl;
  os << "OK=" << S.OK << endl;
 return os;
}


 void IStrahl::initGauss (Vector<complex<double> > &Pol, Gauss g)
{
   double w2,w02,r,zi0,h2,alpha,phi,det,l,R;
 
      Matrix<double> D;
      Vector<double> n,Pkt,Ps,d,hk;

      hk=k;
      d=P-g.F*r0;
      r=sqrt(d[0]*d[0]+d[1]*d[1]);
      w02=g.w0*g.w0*r0*r0;
      zi0=real(k0)*w02/2.0;
      h2=real(k0)*w02/(2*d[2]);
      R=d[2]*(1+h2*h2);

      // Normale zur Phasenfläche
      n[0]=real(k0)*d[0]/R;
      n[1]=real(k0)*d[1]/R;
      n[2]=-1/(zi0*(1+d[2]*d[2]))-real(k0)+real(k0)*r*r/(R*R)*(1-real(k0)*real(k0)*w02*w02/(2*d[2]*d[2]));
      k=n/abs(n);

      // Nächsten Punkt auf der Kugel berechnen
      det=(P*k)*(P*k)-(abs2(P)-r0*r0);
      if (det>0) // Schnittpunkt gefunden
      {
       l=(-P*k-sqrt(det));
       P=P+l*k;
       h2=2*d[2]/(real(k0)*w02);
       w2=w02*(1+h2*h2);
      alpha=abs(ez*k);
      phi=atan(d[2]/zi0)-real(k0)*(d[2]+r*r/(2*R));
      E1=Pol*g.w0*r0/w2*exp(-r*r/w2)*exp(I*phi);      
      D=drehmatrix(k%hk,alpha);
      E1=D*E1;
      E2=E1;
      }
}

 void IStrahl::initEfeldFokus (double sigma2, Vector<double> focuspos,  Vector<complex<double> > Pol)
 { 
   k=(focuspos-P);
   k=k/abs(k);

   double r=(focuspos-P)*k;
   double h=r-abs(focuspos-P);
   
   E1=Pol*exp(-I*k0*h)*exp(-(r*r-h*h)/sigma2);
   E2=E1;
 }

Vector<double> IStrahl::intersectRect (const Vector<double> Pe, const Vector<double> e1, const Vector<double> e2)
{
 Matrix<double> M(e1,e2,-k);
 bool inv;
 cout << "M=" << M << endl;
 Matrix<double>  Mi=invert(M,inv);
 if (inv) 
 {
 Vector<double> Ps=P-Pe;
 Vector<double> L=Mi*Ps;
 cout << "L=" << L << "   k=" << k << "   e1=" << e1 << "   e2="  << e2 << endl;
  if ((L[0]<0.0) || (L[0]>1.0) || (L[1]<0.0) || (L[1]>1.0)) return nanV(""); 
 return P+L[2]*k;
 } 
}
