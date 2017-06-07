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
#include "kegel.h"

Kegel::Kegel()
 : Form()
{
}

Kegel::Kegel (const Form &F) : Form (F)
{
 type = KREISKEGEL;
}

Kegel::Kegel (const Kegel &E)
:Form (E)
{
 r=E.r;
 h=E.h; 
 type = KREISKEGEL;
}


Kegel::~Kegel()
{
}

Kegel::Kegel (
             const Vector<double> &P,
             double r,
             double h,
             complex<double>  n,
             double r0,
             const Matrix<complex<double> > alpha,
             const Vector<double> &Ex,
             const Vector<double> &Ey,
             const Vector<double> &Ez
             ) : Form (P, n, alpha, Ex, Ey, Ez, KREISKEGEL)
{
 this->r0=r0;
 this->h=h;
 this->r=r; 
}



/*!
    \fn Kegel::scale (double sf)
 */
void Kegel::scale (double sf)
{
    h=h/this->sf*sf;
    r=r/this->sf*sf;
    this->sf=sf;
}



/*!
    \fn Kegel::next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside)
 */
bool Kegel::next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
 Vector<double> k,p=Ps-P; 
 double LM,LD,L;
 p=H*p;
 k=H*K;
 LM=schneideMantel (p,k);
 LD=schneideDeckel (p,k);  
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
    \fn Kegel::norm(const Vector<double> &P)
 */
Vector<double> Kegel::norm(const Vector<double> &Ps)
{
 Vector<double> hv,p=H*(Ps-P);
 if (fabs(p[2])<EPS) // Punkt liegt auf Deckel 
 { 
  return -R*ez;   
 }
 else
 { 
   Vector<double> n1,n2;
   n1=p%(h*ez);
   n2=h*ez-p;
   hv=n2%n1;
   hv=hv/abs(hv);
   return R*hv; 
   }
 
 return hv;
}



/*!
    \fn Kegel::isInside(const Vector<double> &Ps)
 */
bool Kegel::isInside(const Vector<double> &Ps)
{
 Vector<double> p=H*(Ps-P);
 return (fabs(p[2])<=r) && (p[0]*p[0]+p[1]*p[1]<r*r);
}



/*!
    \fn Kegel::Volume()
 */
double Kegel::Volume()
{
 return M_PI*r*r*h/3.0;
}

void Kegel::binWrite (ofstream &os)
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

void Kegel::binRead (ifstream &is)
{
 type=ELLIPSOID;
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
    \fn Kegel::setr0(double r0)
 */
void Kegel::setr0(double r0)
{
  r=r/this->r0*r0;
  h=h/this->r0*r0;
  this->r0=r0;
}


/*!
    \fn Kegel::setParms(double r, double h)
 */
void Kegel::setParms(double r, double h)
{
 this->r=r;
 this->h=h;
}


/*!
    \fn Kegel::schneideDeckel(const Vector<double> &P, const Vector<double> &k))
   Annahme: Transformation in das Kegel-Koordinatensystem wurde schon durchgeführt !
            k[2]!=0
 */
double Kegel::schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP)
{
  // erst mal Schnittpunkt mit der Deckelebene berechnen 
   double L=-P[2]/k[2];
    if (L<1E-10*r0) return -1; 
    SP=P+L*k;
   if (SP[0]*SP[0]+SP[1]*SP[1]>r*r) return -1;
   return L;   
}


/*!
    \fn Kegel::schneideMantel (Vector<double> &p, Vector<double> &k)
 */
double Kegel::schneideMantel (Vector<double> &p, Vector<double> &k)
{
  double A,B,C,L1,L2;
  double r2=r*r;
  double h2=h*h;
  A=-r2/h2*k[2]*k[2]+k[0]*k[0]+k[1]*k[1];
  B=2.0*(p[0]*k[0]+p[1]*k[1]-r2/h2*k[2]*p[2]+r2/h*k[2]);
  C=p[0]*p[0]+p[1]*p[1]-r2-r2/h2*p[2]*p[2]+2.0*r2/h*p[2];
  double D=B*B-4.0*A*C;
  double L;
  if (D<0) return -1;  // Es gibt keinen Schnittpunkt
  L1=(-B-sqrt(D))/(2.0*A);
  L2=(-B+sqrt(D))/(2.0*A);
  if (L1>L2) {L=L1; L1=L2; L2=L;} 
  if (L1<EPS) L=L2;
  else L=L1;
  if (L<EPS) L=(-B+sqrt(D))/(2.0*A);  
  Vector<double> t=p+L*k;
  if ((t[2]>h) || (t[2]<0) ) L=-1;    
  return L;
}


/*!
    \fn Kegel::schneideDeckel(Vector<double> &p, Vector<double> &k)
 */
double Kegel::schneideDeckel(Vector<double> &p, Vector<double> &k)
{   
  double L=-p[2]/k[2];
  if (L<EPS)return -1;
   Vector <double> t=p+L*k;  
  if (t[0]*t[0]+t[1]*t[1]>r*r+EPS*EPS) L=-1;
  return L;
}

void Kegel::initQuad() 
{
 
};