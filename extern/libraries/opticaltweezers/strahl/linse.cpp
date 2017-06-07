#include "linse.h"

Linse::Linse()
{
  h=0;
  type = LINSE;
}

Linse::Linse(complex<double> n, Vector<double> P, Vector<double> N, double r1, double r2, double l1, double l2, double h)
{
 this->n=n;
 this->r[0]=r1;
 this->r[1]=r2;
 this->Pl=P;
 this->h=sqrt(r1*r1-l1*l1);
 l[0]=l1;
 l[1]=l2;
 if (this->h>h) this->h=h;
 this->N=N;
 this->P[0]=Pl-l1*N;
 this->P[1]=Pl+l2*N;
 this->r2[0]=r[0]*r[0];
 this->r2[1]=r[1]*r[1];
 this->type = LINSE;
 /*cout << "P0=" << this->P[0] << "     P1=" << this->P[1] << endl;
 cout << "l0=" << this->l[0] << "    l1=" << this->l[1] << endl;
 cout << "N=" << N << endl;*/
}

void Linse::scale(double sf)
{
  this->sf=sf;
}


#define sign(x) ((x<0) ? 1 : 0) 

bool Linse::next(const Vector<double> &p,const  Vector<double> &k, 
                                 Vector <double> &pout, const  int inside)
{
 int i;
 double l0,l1,l2,B;
 Vector<double> Psi;
 bool innen=(sgn(l[0])*(abs(p-P[0])-r[0])<1E-5) && (sgn(l[1])*(abs(p-P[1])-r[1])<1E-5);
 innen=innen && (abs(p-Pl)<r[0]) && (abs(p-Pl)<r[1]);

 if (innen  ) // Innen
 {
  i=sign(N*k);   
 // cout << "innen" << endl;
  }
 else // außen
 {
  i=sign(-N*k);
 // cout << "außen" << endl;
  }
 /*cout << "|p-P[0]|-r0=" <<sgn(l[0])*( abs(p-P[0])-r[0] )<< "    |p-P[1]|-r1=" << sgn(l[1])*(abs(p-P[1])-r[1]) << endl;
   cout << "i=" << i << "      P=" << P[i] << endl;   */
  Psi=p-P[i];
  B=k*Psi;
  
  
  l0=B*B-Psi*Psi+r2[i];
 
  if (l[i]>0) // Seite ist konvex
  {
  
  if (l0<0) return false; 
  l2=-B+sqrt(l0);
  l1=-B-sqrt(l0);
//   cout << "l1=" << l1 << "    l2=" << l2 << endl; 
  if (l1>0.0) { pout=p+l1*k; return true; }
  else if (innen)
   { 
    if (l2<0.0) return false;  
   pout=p+l2*k; 
   return true;
   }
   return false; 
  }
  else // Seite ist konkav
  {
   l1=-B-sqrt(l0); 
   l2=-B+sqrt(l0);
 //  cout << "l1=" << l1 << "    l2=" << l2 << endl; 
   if ((l1>0) && innen) {pout=p+l1*k; return true; }
   
   if (l2<0.0) return false;
   pout=p+l2*k; 
   return true;
  }
  
  return false;
}

Vector<double> Linse::norm(const Vector<double> &p)
{
  Vector<double> h=p-P[0];
  if (fabs(h*h-r2[0])>1E-10) h=p-P[1]; 
  return h/abs(h);                      // Muss geändert werden
}

bool Linse::isInside(const Vector<double> &p)
{
  Vector<double> h=p-Pl;
  double h2=h*h;
  return (h2-r2[0]<0) && (h2-r2[1]<0);
}

void Linse::setr0(double r0)
{
	Pl=Pl/this->r0*r0;
	P[0]=P[0]/this->r0*r0;
    P[1]=P[1]/this->r0*r0;
	r[0]=r[0]/this->r0*r0;
    r[1]=r[1]/this->r0*r0;
	r2[0]=r2[0]/this->r0*r0;
    r2[1]=r2[1]/this->r0*r0;
	l[0]=l[0]/this->r0*r0;
	l[1]=l[1]/this->r0*r0;
	h=h/this->r0*r0;
	this->r0=r0;
}

void Linse::setP (Vector<double> r) 
{ 
 Pl=r; 
 P[0]=Pl-l[0]*N; 
 P[1]=Pl+l[1]*N;  
 initQuad();
} //



ostream &operator << (ostream &os, Linse L)
{
  os << "P0=" << L.P[0] << "   P1=" << L.P[1] << "    Pl=" << L.Pl << endl;
  os << "n=" << L.N << "   r0=" << L.r[0]  << "   r1=" << L.r[1] << "   h=" << L.h << endl;
  return os; 
}