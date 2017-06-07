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
#include "zylinder.h"
#include "intersect.h"
#include <algorithm>

Zylinder::Zylinder()
 : Form()
{
	type=ZYLINDER;
}

Zylinder::Zylinder (const Form &F) : Form (F)
{
 type = ZYLINDER;
}

Zylinder::Zylinder (const Zylinder &E)
:Form (E)
{
 r=E.r;
 h=E.h; 
 type = ZYLINDER;
}



Zylinder::Zylinder (
             const Vector<double> &P,
             double r,
             double h,
             complex<double>  n,
             double r0,
             const Matrix<complex<double> > alpha,
             const Vector<double> &Ex,
             const Vector<double> &Ey,
             const Vector<double> &Ez
             ) : Form (P, n, alpha, Ex, Ey, Ez,ZYLINDER)
{
 this->r0=r0;
 this->h=h;
 this->r=r; 
}



/*!
    \fn Zylinder::scale (double sf)
 */
void Zylinder::scale (double sf)
{
    h=h/this->sf*sf;
    r=r/this->sf*sf;
    this->sf=sf;
}



/*!
    \fn Zylinder::next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside)
 */
bool Zylinder::next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
 Vector<double> Sp,k,p=Ps-P;  
 double LM,LD,L;
 p=H*p;
 k=H*K;
 LM=schneideMantel (p,k);
 LD=schneideDeckel (p,k,Sp);  
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




/*!
    \fn Zylinder::norm(const Vector<double> &P)
 */
Vector<double> Zylinder::norm(const Vector<double> &Ps)
{
 Vector<double> hv,p=H*(Ps-P);
 if (p[2]<EPS) hv=-R*ez;
 else
 {
   if (p[2]>h-EPS) hv=R*ez;  
   else hv=R*Vector<double> (p[0],p[1],0); 
 }    
  hv=hv/abs(hv);
 return hv;
}



/*!
    \fn Zylinder::isInside(const Vector<double> &Ps)
 */
bool Zylinder::isInside(const Vector<double> &Ps)
{
 Vector<double> p=H*(Ps-P);
 return (fabs(p[2])<=r) && (p[0]*p[0]+p[1]*p[1]<r*r);
}



/*!
    \fn Zylinder::Volume()
 */
double Zylinder::Volume()
{
 return M_PI*r*r*h;
}

void Zylinder::binWrite (ofstream &os)
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
 os.write ((char *) &r,(char) sizeof (r));
 os.write ((char *) &h,(char) sizeof (h));
}

void Zylinder::binRead (ifstream &is)
{
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
 is.read ((char *) &r,(char) sizeof (r));
 is.read ((char *) &h,(char) sizeof (h)); 
}


/*!
    \fn Zylinder::setr0(double r0)
 */
void Zylinder::setr0(double r0)
{
  r=r/this->r0*r0;
  h=h/this->r0*r0;
  P=P/this->r0*r0;
  this->r0=r0;  
}


/*!
    \fn Zylinder::setParms(double r, double h)
 */
void Zylinder::setParms(double r, double h)
{
 this->r=r;
 this->h=h;
}


/*!
    \fn Zylinder::schneideDeckel(const Vector<double> &P, const Vector<double> &k))
    findet Schnittpunkt mit einem der Deckel
   Annahme: Transformation in das Zylinder-Koordinatensystem wurde schon durchgeführt !
            k[2]!=0 
     Rückgabewert: Lauflänge vom Momentanen Punkt P zum Schnittpunkt in Strahlrichtung  
 */
double Zylinder::schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP)
{
  // erst mal Schnittpunkt mit der Deckelebene berechnen 
   double L1=-P[2]/k[2];
   double L2=(h-P[2])/k[2];
   double L;
/*    if (L1<r0*1E-10) L=L2; 
  else 
    { 
     if (L1<L2) L=L1;
     else L=L2;
     }
  */ 
   
   if ((L1<L2) && (L1>EPS)) L=L1;
   else 
	   if (L2>EPS) L=L2;
	   else L=L1; 
   if (L<r0*1E-20) return -1;
   SP=P+L*k;
   if (SP[0]*SP[0]+SP[1]*SP[1]>r*r) return -1;
   return L;   
}


/*!
    \fn Zylinder::schneideMantel (Vector<double> &p, Vector<double> &k)
 */
double Zylinder::schneideMantel (Vector<double> &p, Vector<double> &k)
{
  double A,B,C;
  A=k[0]*k[0]+k[1]*k[1];
  B=2.0*(p[0]*k[0]+p[1]*k[1]);
  C=p[0]*p[0]+p[1]*p[1]-r*r;
  double D=B*B-4.0*A*C;
  double L;
  if (D<0) return -1;  // Es gibt keinen Schnittpunkt
  L=(-B-sqrt(D))/(2.0*A);
  if (L<EPS) L=(-B+sqrt(D))/(2.0*A);  
  Vector<double> t=p+L*k;
  if ((t[2]>h-EPS) || (t[2]<EPS)) L=-1;    
  return L;
}


/*!
    \fn Zylinder::schneideDeckel(Vector<double> &p, Vector<double> &k)
 */
double Zylinder::schneideDeckel(Vector<double> &p, Vector<double> &k)
{   
	 cout << "HÄH?" << endl;
  double L=(-h/2.0-p[2])/k[2];
  if (L<EPS) L=L=(h/2.0-p[2])/k[2];
   Vector <double> t=p+L*k;  
  if (t[0]*t[0]+t[1]*t[1]>r*r+EPS*EPS) L=-1;
  return L;
}

void Zylinder::initQuad()
{
   Vector<double> hv[8];

   for (int k=0; k<=1; k++)
   for (int l=0; l<=1; l++)
     for (int m=0; m<=1; m++)
     {
        k==0 ? hv[k*4+l*2+m][0]=-r : hv[k*4+l*2+m][0]=r;
        l==0 ? hv[k*4+l*2+m][1]=-r : hv[k*4+l*2+m][1]=r;
        m==0 ? hv[k*4+l*2+m][2]=0 : hv[k*4+l*2+m][2]=h;
     }

 pul=zero;
 por=zero;
 for (int i=0; i<8; i++)
 {
  hv[i]=R*hv[i];
  if (hv[i][0]<pul[0]) pul[0]=hv[i][0];
  if (hv[i][1]<pul[1]) pul[1]=hv[i][1];
  if (hv[i][2]<pul[2]) pul[2]=hv[i][2];
  if (hv[i][0]>por[0]) por[0]=hv[i][0];
  if (hv[i][1]>por[1]) por[1]=hv[i][1];
  if (hv[i][2]>por[2]) por[2]=hv[i][2];
 }
 pul=pul+P;
 por=por+P;
 cout << "Zylinder : pul=" << pul << "    por=" << por << endl;
}


ZylinderHexagonal::ZylinderHexagonal(double a, double h)
{
 this->a=a;
 this->h=h;
 type=ZYLINDER_HEXAGONAL;
 init();
}

void ZylinderHexagonal::init()
{
 double phi;
 double alpha=1.0/3.0*M_PI; // entspricht 60°
 for (int i=0; i<6; i++)
 {
  phi=(i+2)*alpha;
  e1[i]=Vector<double> (a*cos((i+2)*alpha),a*sin((i+2)*alpha),0);
  Pe[i]=Vector<double> (a*cos(i*alpha),a*sin(i*alpha),0);
  ne[i]=Vector<double> (cos((i+0.5)*alpha),sin((i+0.5)*alpha),0);
 }
}
 
ostream& operator << (ostream &os, ZylinderHexagonal Z)
{
  int i;
  os << "type=" << Z.type << endl;
  os << "Pe:" << endl;
  for (i=0; i<6; i++)
    os << Z.Pe[i] << endl;
  os << "e1:" << endl;
  for (i=0; i<6; i++)
    os << Z.e1[i] << endl;
  os << "ne:" << endl;
  for (i=0; i<6; i++)
    os << Z.ne[i] << endl;
 return os;  
}

double ZylinderHexagonal::schneideMantel (Vector<double> &p, Vector<double> &k)
{
  int i,n=0;
  Vector<double> P[2],L;
  double l[2];
  l[0]=nan("");
  l[1]=nan("");
//  cout << "p=" << p << "  k=" << k << endl;
 for (i=0; (i<6) && (n<2); i++)
  {
     P[n]=intersectLineRect(Pe[i],e1[i],ez*h,p,k,L);
 //    cout << "i=" << i << "  L=" << L << endl; 
     if ((!isnan(P[n])) && (L[2]>0.0)) 
          {
           l[n]=L[2];
           n++;
           }  
  }
//   cout << "l0=" << l[0] << "   l1=" << l[1] << endl;
  return std::min(l[0],l[1]);
}

