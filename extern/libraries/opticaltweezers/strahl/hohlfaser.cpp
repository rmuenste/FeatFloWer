#include "hohlfaser.h"

Hohlfaser::Hohlfaser(void)
{ 
	 type=HOHLFASER;
}


Hohlfaser::Hohlfaser (const Hohlfaser &E)
:Zylinder ((Zylinder) E)
{
 ra=E.ra;
 ri=E.ri;
 h=E.h; 
 type = HOHLFASER;
}

   
Hohlfaser::Hohlfaser (
             const Vector<double> &P,
             double ri,
			 double ra,
             double h,
             complex<double>  n,
             double r0,
             const Matrix<complex<double> > alpha,
             const Vector<double> &Ex,
             const Vector<double> &Ey,
             const Vector<double> &Ez
             ) : Zylinder (P, ra, h, n, r0, alpha, Ex, Ey, Ez)
{
 this->r0=r0;
 this->h=h;
 this->ra=ra;
 this->ri=ri;
 r=ra;
 type = HOHLFASER;
}


Hohlfaser::~Hohlfaser(void)
{
}

double Hohlfaser::schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP)
{
    double L;
	L=Zylinder::schneideDeckel(P,k,SP);
	if (sqrt(SP[0]*SP[0]+SP[1]*SP[1])<ri) {SP=P; return -1; }
	else return L;
}


double Hohlfaser::schneideMantel (Vector<double> &p, Vector<double> &k)
{
  double A,B,C,D;
  Vector<double> t;
  double L,L1,L2;

  // erst mal die innere Mantelflaeche
  A=k[0]*k[0]+k[1]*k[1];
  if (A==0) return -1;
  B=2.0*(p[0]*k[0]+p[1]*k[1]);
  C=p[0]*p[0]+p[1]*p[1]-ri*ri;
  D=B*B-4.0*A*C;
  
  if (D<0) L1=-1;  // Es gibt keinen Schnittpunkt
  else
  {	 
  L1=(-B-sqrt(D))/(2.0*A);
  if (L1<EPS) L1=(-B+sqrt(D))/(2.0*A);  
  t=p+L1*k;
  if ((t[2]>h-EPS) || (t[2]<EPS)) L1=-1;    
  }


   // jetzt die aeussere Mantelflaeche
  /*A=k[0]*k[0]+k[1]*k[1];
  B=2.0*(p[0]*k[0]+p[1]*k[1]);*/ // <- nicht doppelt zu berechnen
  C=p[0]*p[0]+p[1]*p[1]-ra*ra;
  D=B*B-4.0*A*C;
  
  if (D<0) L2=-1;  // Es gibt keinen Schnittpunkt
  else
  {	 
  L2=(-B-sqrt(D))/(2.0*A);
  if (L2<EPS) L2=(-B+sqrt(D))/(2.0*A);  
  t=p+L1*k;
  if ((t[2]>h-EPS) || (t[2]<EPS)) L2=-1;    
  }

  if ((L1<L2) && (L1>0)) L=L1;
  else 
     if (L2<0) L=-1;
  	 else L=L2;
  return L; 	 
}

bool Hohlfaser::next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
 Vector<double> Sp,k,p=Ps-P;  
 double LM,LD,L;
 p=H*p;
 k=H*K;
 LM=schneideMantel (p,k);
 cout << "LM=" << LM << endl;
 LD=schneideDeckel (p,k,Sp); 
 cout << "LD=" << LD << endl;
 if (LM>LD) 
 {
  if (LD>EPS) L=LD; 
  else L=LM;
 }
 else 
 if (LM>EPS) L=LM;
 else L=LD;

 if (L<EPS) {pout=Ps; return false;}
 pout=Ps+L*K;
 return true; 
}

Vector<double> Hohlfaser::norm(const Vector<double> &Ps)
{
 Vector<double> hhv,hv,p=H*(Ps-P);
 double h,hr;
 h=0;
 if (p[2]<EPS) hv=-R*ez;
 else
 {
   if (p[2]>h-EPS) hv=R*ez;  
   else
   {
	   hhv=Vector<double> (p[0],p[1],0); 
	   hr=abs(hhv);
	   if (hr<ra-EPS) hv=-hhv;
	   else hv=hhv;
   }
 }    
  hv=hv/abs(hv);
 return hv;
}
