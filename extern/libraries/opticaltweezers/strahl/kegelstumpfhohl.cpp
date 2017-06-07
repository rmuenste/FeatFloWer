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
#include "kegelstumpfhohl.h"

KegelstumpfHohl::KegelstumpfHohl(const Vector< double >& P, 
    double rimin, double rimax, double ramin, double ramax, 
    double h, complex< double > n, 
    double r0, 
    const Matrix< complex < double > > alpha, 
    const Vector< double >& Ex, 
    const Vector< double >& Ey, 
    const Vector< double >& Ez)
 : Form(P, n, alpha, Ex, Ey, Ez, KEGELSTUMPF_HOHL)
{
/* if (rimin>rimax) {this->rimin=rimax; this->rimin=rimin; }
 else { this->rimin=rimin; this->rimax=rimax; };*/

 this->rimin=rimin; this->rimax=rimax; 
 /*if (ramin>ramax) {this->ramin=ramax; this->ramin=ramin; }
 else { this->ramin=ramin; this->ramax=ramax; }; */
 this->ramin=ramin; this->ramax=ramax; 
 this->h=h;
 hsa=this->ramax/(this->ramax-this->ramin)*this->h;
 hsi=this->rimax/fabs((this->rimax-this->rimin))*this->h;
 this->r0=r0;
 type=KEGELSTUMPF_HOHL; 
}

KegelstumpfHohl::KegelstumpfHohl (const KegelstumpfHohl &E) : Form (E)
{
 this->ramin=E.ramin;
 this->ramax=E.ramax;
 this->rimin=E.rimin;
 this->rimax=E.rimax;
 this->h=E.h;
 this->hsa=E.hsa;
 this->hsi=E.hsi;
 type=KEGELSTUMPF_HOHL;
} 

KegelstumpfHohl::KegelstumpfHohl (const Form &E) : Form (E)
{
 type=KEGELSTUMPF_HOHL;
} 

double KegelstumpfHohl::Volume()
{
 double V1i,V2i,V1a,V2a;
 double Vi,Va;
 double l=hsi-h; 
 V1i=M_PI*rimax*rimax*hsi/3.0;
 V2i=M_PI*rimin*rimin*l/3.0;
 Vi=V1i-V2i;

 l=hsa-h;
 V1a=M_PI*ramax*ramax*hsa/3.0;
 V2a=M_PI*ramin*ramin*l/3.0;
 Va=V1a-V2a;  
 return Va-Vi; 
}

double KegelstumpfHohl::schneideDeckel(int i, Vector<double> &p, Vector<double> &k)
{ 
  double L;
  if (k[2]==0) return -1;
  if (i==0) // unterer Deckel
  {   
   L=-p[2]/k[2];
   if (L<EPS) return -1;
   Vector <double> t=p+L*k;  
   if ((t[0]*t[0]+t[1]*t[1]<(rimax*rimax-EPS*EPS)) || (t[0]*t[0]+t[1]*t[1]>(ramax*ramax+EPS*EPS)) ) L=-1;
  }
  else  // oberer Deckel
  {
   L=(h-p[2])/k[2];
   if (L<EPS) return -1;
   Vector <double> t=p+L*k;  
   if ((t[0]*t[0]+t[1]*t[1]<(rimin*rimin-EPS*EPS)) || (t[0]*t[0]+t[1]*t[1]>(ramin*ramin+EPS*EPS)) ) L=-1;
  }
  return L;
}

double KegelstumpfHohl::schneideMantel (int i, Vector<double> &p, Vector<double> &k)
{
  Vector<double> t;
  double A,B,C,L1,L2;
  double r2;
  double D,L,h2;
  if (i==0) // innerer Mantel
  { 
   r2=rimax*rimax;
   h2=hsi*hsi;
   A=-r2/h2*k[2]*k[2]+k[0]*k[0]+k[1]*k[1];
   B=2.0*(p[0]*k[0]+p[1]*k[1]-r2/h2*k[2]*p[2]+r2/hsi*k[2]);
   C=p[0]*p[0]+p[1]*p[1]-r2-r2/h2*p[2]*p[2]+2.0*r2/hsi*p[2];
   D=B*B-4.0*A*C;   
   if (D<0) return -1;  // Es gibt keinen Schnittpunkt
   L1=(-B-sqrt(D))/(2.0*A);
   L2=(-B+sqrt(D))/(2.0*A);  
   if (L1>L2) {L=L1; L1=L2; L2=L;} 
   if (L1<EPS) L=L2;
   else L=L1;
   if (L<EPS) L=(-B+sqrt(D))/(2.0*A);
   t=p+L*k;
   if ((t[2]>h) || (t[2]<0)) L=-1;
  }
  else
  { 
   r2=ramax*ramax;
   h2=hsa*hsa;
   A=-r2/h2*k[2]*k[2]+k[0]*k[0]+k[1]*k[1];
   B=2.0*(p[0]*k[0]+p[1]*k[1]-r2/h2*k[2]*p[2]+r2/hsa*k[2]);
   C=p[0]*p[0]+p[1]*p[1]-r2-r2/h2*p[2]*p[2]+2.0*r2/hsa*p[2];
   D=B*B-4.0*A*C;
   L;
   if (D<0) return -1;  // Es gibt keinen Schnittpunkt
   L1=(-B-sqrt(D))/(2.0*A);
   L2=(-B+sqrt(D))/(2.0*A);  
   if (L1>L2) {L=L1; L1=L2; L2=L;} 
   if (L1<EPS) L=L2;
   else L=L1;
   if (L<EPS) L=(-B+sqrt(D))/(2.0*A);
   t=p+L*k;
   if ((t[2]>h) || (t[2]<0)) L=-1;
  }  
  return L;
}

bool KegelstumpfHohl::next(const Vector< double >& Ps, const Vector< double >&K, Vector< double >& pout, const int inside)
{
 Vector<double> k,p=Ps-P; 
 double LM,LH[4]={-1,-1,-1,-1};
 double L;
 p=H*p;
 k=H*K;
 LH[0]=schneideMantel (0,p,k);  
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


 LM=schneideMantel (1,p,k);
 if (LM>LH[2]) {LH[3]=LM;}
 else 
 { 
  if (LM>LH[1]) { LH[3]=LH[2];  LH[2]=LM; }
  else 
    { 
      if (LM>LH[0]) { LH[3]=LH[2]; LH[2]=LH[1]; LH[1]=LM; }
      else { LH[3]=LH[2]; LH[2]=LH[1]; LH[1]=LH[0]; LH[0]=LM; }
    } 
 }

 L=-1;
 for (int i=0; (i<4) && (L<EPS); i++)
   L=LH[i]; 
 if (L<EPS) {pout=Ps; return false;}
 pout=Ps+L*K;
 return true; 
}

Vector<double> KegelstumpfHohl::norm(const Vector<double> &Ps)
{
 Vector<double> hv,p=H*(Ps-P);
Vector<double> n1,n2;
 double r;
 if (fabs(p[2])<EPS) // Punkt liegt auf Deckel 
 { 
  return -R*ez;   
 }

 if (fabs(p[2])>=h-EPS) // Punkt liegt auf Deckel 
 { 
  return R*ez;   
 }
 r=sqrt(p[0]*p[0]+p[1]*p[1]);
 if (fabs(p[0]+rimax/hsa*r-hsi)>EPS) // Punkt liegt auf dem ‰uﬂeren Mantel
 {  
   n1=p%(hsa*ez);
   n2=hsa*ez-p;
 }
 else 
 {
   n1=p%(hsi*ez);
   n2=hsi*ez-p;
 } 
   hv=n2%n1;
   hv=hv/abs(hv);
   return R*hv; 
}


KegelstumpfHohl::~KegelstumpfHohl()
{
}

void KegelstumpfHohl::setParms (double rimin, double rimax, double ramin, double ramax, double h)
{
 /*if (rimin<rimax) { this->rimin=rimin; this->rimax=rimax; }
 else { this->rimax=rimin; this->rimax=rimin; }*/

 this->rimin=rimin; this->rimax=rimax; 
 if (ramin<ramax) { this->ramin=ramin; this->ramax=ramax; }
 else { this->ramax=ramin; this->ramax=ramin; }

 this->h=h;
 hsa=h*this->ramax/(this->ramax-this->ramin);
 hsi=h*this->rimax/(this->rimax-this->rimin);
}

void KegelstumpfHohl::binWrite (ofstream &os)
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
 os.write ((char *) &rimin,(char) sizeof (rimin));
 os.write ((char *) &rimax,(char) sizeof (rimax));
 os.write ((char *) &ramin,(char) sizeof (ramin));
 os.write ((char *) &ramax,(char) sizeof (ramax));
 os.write ((char *) &h,(char) sizeof (h));
}

void KegelstumpfHohl::binRead (ifstream &is)
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
 is.read ((char *) &rimin,(char) sizeof (rimin));
 is.read ((char *) &rimax,(char) sizeof (rimax));
 is.read ((char *) &ramin,(char) sizeof (ramin));
 is.read ((char *) &ramax,(char) sizeof (ramax));
 is.read ((char *) &h,(char) sizeof (h)); 
}

void KegelstumpfHohl::scale (double sf)
{
    h=h/this->sf*sf;
    hsa=hsa/this->sf*sf;
    hsi=hsi/this->sf*sf;
    rimin=rimin/this->sf*sf;
    rimax=rimax/this->sf*sf;
    ramin=ramin/this->sf*sf;
    ramax=ramax/this->sf*sf;
    this->sf=sf;
}

void KegelstumpfHohl::setr0(double r0)
{
  rimin=rimin/this->r0*r0;
  rimax=rimax/this->r0*r0;
  ramin=ramin/this->r0*r0;
  ramax=ramax/this->r0*r0; 
  h=h/this->r0*r0;
  hsa=hsa/this->r0*r0;
  hsi=hsi/this->r0*r0; 
  this->r0=r0;
}
