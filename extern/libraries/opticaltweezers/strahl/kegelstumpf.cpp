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
#include "kegelstumpf.h"

Kegelstumpf::Kegelstumpf(const Vector< double >& P, double rmin, double rmax, double h, complex< double > n, double r0, const Matrix< complex < double > > alpha, const Vector< double >& Ex, const Vector< double >& Ey, const Vector< double >& Ez) : Form (P, n, alpha, Ex, Ey, Ez, KEGELSTUMPF)
{ 
 if (rmin>rmax) {this->rmax=rmin; this->rmin=rmax; }
 else
 {
 this->rmin=rmin;
 this->rmax=rmax;
 }  
 this->h=h;  
 hs=this->rmax/(this->rmax-this->rmin)*this->h;
 this->r0=r0;
 cout << "rmin=" << this->rmin << "   rmax=" << this->rmax << "   hs=" << this->hs << endl; 
 type=KEGELSTUMPF;
}


Kegelstumpf::Kegelstumpf(const Form& F): Form(F)
{
 type=KEGELSTUMPF; 
}



Kegelstumpf::Kegelstumpf (const Kegelstumpf &E)
:Form (E)
{
 r0=E.r0;
 rmin=E.rmin;
 rmax=E.rmax;
 h=E.h;
 hs=E.hs;
 type = KEGELSTUMPF;
}




Kegelstumpf::~Kegelstumpf()
{
}


bool Kegelstumpf::isInside(const Vector< double >& Ps)
{
    return true; // MUSS NOCH geändert werden
}

bool Kegelstumpf::next(const Vector< double >& Ps, const Vector< double >&K, Vector< double >& pout, const int inside)
{
 Vector<double> k,p=Ps-P; 
 double LM,LH[3]={-1,-1,-1};
 double L;
 p=H*p;
 k=H*K;
 LH[0]=schneideMantel (p,k);  
 LM=schneideDeckel (0,p,k);
 if (LM<LH[0]) {LH[1]=LH[0]; LH[0]=LM; }
 else LH[1]=LM;  
 LM=schneideDeckel (1,p,k);
 if (LM>LH[1]) { LH[2]=LM; }
 else
 {
  if (LM>LH[0]) { LH[2]=LH[1]; LH[1]=LM; }
  else  { LH[2]=LH[1]; LH[1]=LH[0]; LH[0]=LM; }
 }
 
 L=-1;
 for (int i=0; (i<3) && (L<EPS); i++)
   L=LH[i]; 
 if (L<EPS) {pout=Ps; return false;}
 pout=Ps+L*K;
 return true; 
}

double Kegelstumpf::schneideDeckel(int i, Vector<double> &p, Vector<double> &k)
{ 
  double L;
  if (k[2]==0) return -1;
  if (i==0) // unterer Deckel
  {   
   L=-p[2]/k[2];
   if (L<EPS) return -1;
   Vector <double> t=p+L*k;  
   if (t[0]*t[0]+t[1]*t[1]>(rmax*rmax+EPS*EPS)) L=-1;
  }
  else  // oberer Deckel
  {
   L=(h-p[2])/k[2];
   if (L<EPS) return -1;
   Vector <double> t=p+L*k;  
   if (t[0]*t[0]+t[1]*t[1]>(rmin*rmin+EPS*EPS)) L=-1;
  }
  return L;
}

Vector<double> Kegelstumpf::norm(const Vector<double> &Ps)
{
 Vector<double> hv,p=H*(Ps-P);
 if (fabs(p[2])<EPS) // Punkt liegt auf Deckel 
 { 
  return -R*ez;   
 }
 if (fabs(p[2])>=h-EPS) // Punkt liegt auf Deckel 
 { 
  return R*ez;   
 }

 Vector<double> n1,n2;
   n1=p%(hs*ez);
   n2=hs*ez-p;
   hv=n2%n1;
   hv=hv/abs(hv);
   return R*hv; 
}

 
double Kegelstumpf::schneideMantel (Vector<double> &p, Vector<double> &k)
{
  Vector<double> t;
  double A,B,C,L1,L2;
  double r2=rmax*rmax;
  double h2=hs*hs;
  A=-r2/h2*k[2]*k[2]+k[0]*k[0]+k[1]*k[1];
  B=2.0*(p[0]*k[0]+p[1]*k[1]-r2/h2*k[2]*p[2]+r2/hs*k[2]);
  C=p[0]*p[0]+p[1]*p[1]-r2-r2/h2*p[2]*p[2]+2.0*r2/hs*p[2];
  double D=B*B-4.0*A*C;
  double L;
  if (D<0) return -1;  // Es gibt keinen Schnittpunkt
  L1=(-B-sqrt(D))/(2.0*A);
  L2=(-B+sqrt(D))/(2.0*A);  
  if (L1>L2) {L=L1; L1=L2; L2=L;} 
  if (L1<EPS) L=L2;
  else L=L1;
  if (L<EPS) L=(-B+sqrt(D))/(2.0*A);
  t=p+L*k;
  if ((t[2]>h) || (t[2]<0)) L=-1;  
  return L;
}


void Kegelstumpf::binWrite (ofstream &os)
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
 os.write ((char *) &rmin,(char) sizeof (rmin));
 os.write ((char *) &rmax,(char) sizeof (rmax));
 os.write ((char *) &h,(char) sizeof (h));
}

void Kegelstumpf::binRead (ifstream &is)
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
 is.read ((char *) &rmin,(char) sizeof (rmin));
 is.read ((char *) &rmax,(char) sizeof (rmax));
 is.read ((char *) &h,(char) sizeof (h)); 
}


void Kegelstumpf::setr0(double r0)
{
  rmin=rmin/this->r0*r0;
  rmax=rmax/this->r0*r0; 
  h=h/this->r0*r0;
  hs=hs/this->r0*r0;
  this->r0=r0;
}

void Kegelstumpf::scale (double sf)
{
    h=h/this->sf*sf;
    hs=hs/this->sf*sf;
    rmin=rmin/this->sf*sf;
    rmax=rmax/this->sf*sf;
    this->sf=sf;
}

double Kegelstumpf::Volume()
{
 double V1,V2;
 double l=hs-h; 
 V1=M_PI*rmax*rmax*hs/3.0;
 V2=M_PI*rmin*rmin*l/3.0;
 return V1-V2; 
}

void Kegelstumpf::initQuad() 
{
 
};