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
#include "superellipsoid.h"

#ifndef min(x,y)
#define min(x,y) x<y ? x : y
#endif


Superellipsoid_D::Superellipsoid_D()
 : surface()
{
 FName="Superellipsoid_D"; 
 type=SUPERELLIPSOID_D;   
 r0=1.0;
 sf=1.0;
 e1=0;
 e2=0;
 a1=0;
 a2=0;
 a3=0;
 ntheta=0;
 nphi=0;
}

Superellipsoid_D::Superellipsoid_D (const Superellipsoid_D &Su):surface((surface)Su)
{
  FName="Superellipsoid_D";
  anzp=Su.anzp;
  sf=Su.sf;
  r0=Su.r0;
  sf=1.0;  
  P = zero;
  currentnorm=Vector<double>(0,0,0);
  S=new dreieck[anzp];
  P=Su.P;
  for (int i=0;i<anzp;i++)
   S[i]=Su.S[i];
  type=SUPERELLIPSOID_D;  
  a1=Su.a1;
  a2=Su.a2;
  a3=Su.a3;
  e1=Su.e1;
  e2=Su.e2;
  nphi=Su.nphi;
  ntheta=Su.ntheta;
}



Superellipsoid_D::~Superellipsoid_D()
{
}

Vector<double> Superellipsoid_D::super (double theta, double phi,double a1, double a2, double a3, double e1, double e2)
{
 Vector<double> erg;
 double ct=cos(theta);
 double st=sin(theta);
 double sp=sin(phi);
 double cp=cos(phi);

 erg[0]=a1*sign(st)*pow(fabs(st),e1)*sign(cp)*pow(fabs(cp),e2);
 erg[1]=a2*sign(st)*pow(fabs(st),e1)*sign(sp)*pow(fabs(sp),e2);
 erg[2]=a3*sign(ct)*pow(fabs(ct),e1);

 return erg;
}

void Superellipsoid_D::generate (int ntheta, int nphi, double a1, double a2, double a3, double e1, double e2)
{
 this->ntheta=ntheta;
 this->nphi=nphi;
 this->a1=a1;
 this->a2=a2;
 this->a3=a3;
 this->e1=e1;
 this->e2=e2;
 Vector<double> P[4];
 type=SUPERELLIPSOID_D;
 sf=1.0;
 this->clearS();
//  if (ntheta%2!=0) ntheta++; 

 this->n=n;
 double dtheta=M_PI/(double)(ntheta);
 double dphi=2.0*M_PI/(double)(nphi);

 double theta,phi;
 int count=0;
 if ((ntheta>0) && (nphi>0)) 
 {
 anzp=2*(ntheta-2)*nphi + 2*nphi;
 S=new dreieck[anzp];
 
 /* ------- Mittelteil (ohne Deckel) ---------------- */
 for (int itheta=1; itheta<=ntheta-2; itheta++) 
 {
  theta=itheta*dtheta;
  for (int iphi=0; iphi<=nphi-1; iphi++)
  {
   phi=iphi*dphi;
   P[0]=super (theta,phi,a1,a2,a3,e1,e2);
   P[1]=super (theta+dtheta,phi,a1,a2,a3,e1,e2);
   P[2]=super (theta+dtheta,phi+dphi,a1,a2,a3,e1,e2);
   P[3]=super (theta,phi+dphi,a1,a2,a3,e1,e2);

  S[count]=dreieck(P[0]*r0,P[1]*r0,P[3]*r0);
  count++;
  S[count]=dreieck(P[1]*r0,P[3]*r0,P[2]*r0);
  count++;
  }
 }

/* ------------ Oberer Deckel ---------------- */
 for (int iphi=0; iphi<=nphi-1; iphi++)
  {
   phi=iphi*dphi;
   P[0]=super (0,phi,a1,a2,a3,e1,e2);
   P[1]=super (dtheta,phi,a1,a2,a3,e1,e2);
   P[2]=super (dtheta,phi+dphi,a1,a2,a3,e1,e2);
   S[count]=dreieck(P[0]*r0,P[1]*r0,P[2]*r0);
   count++;
  }


/* ----------- Unterer Deckel --------------------- */
  for (int iphi=0; iphi<=nphi-1; iphi++)
  {
   phi=iphi*dphi;
   P[0]=super (M_PI,phi,a1,a2,a3,e1,e2);
   P[1]=super (M_PI-dtheta,phi,a1,a2,a3,e1,e2);
   P[2]=super (M_PI-dtheta,phi+dphi,a1,a2,a3,e1,e2);
   S[count]=dreieck(P[0]*r0,P[1]*r0,P[2]*r0);
   count++;
  }
//   exportSRF ("test.srf");
//   cout << "volumen=" << calcVolume() << endl;
 }
}

void Superellipsoid_D::binWrite (ofstream &os)
{
 size_t strl;
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
 os.write ((char *) &anzp, (char) sizeof(anzp));
 os.write ((char *) &a1, (char) sizeof(a1));
 os.write ((char *) &a2, (char) sizeof(a2));
 os.write ((char *) &a3, (char) sizeof(a3));
 os.write ((char *) &e1, (char) sizeof(e1));
 os.write ((char *) &e2, (char) sizeof(e2));
 os.write ((char *) &ntheta, (char) sizeof(ntheta));
 os.write ((char *) &nphi, (char) sizeof(nphi));

 for (int i=0; i<anzp; i++)
  S[i].binWrite(os);  
}

void Superellipsoid_D::binRead (ifstream &is)
{
 type=SUPERELLIPSOID_D;
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
 clearS();
 is.read ((char *) &anzp,(char) sizeof (anzp));
 is.read ((char *) &a1, (char) sizeof(a1));
 is.read ((char *) &a2, (char) sizeof(a2));
 is.read ((char *) &a3, (char) sizeof(a3));
 is.read ((char *) &e1, (char) sizeof(e1));
 is.read ((char *) &e2, (char) sizeof(e2));
 is.read ((char *) &ntheta, (char) sizeof(ntheta));
 is.read ((char *) &nphi, (char) sizeof(nphi)); 
 S=new dreieck[anzp];
 for (int i=0; i<anzp; i++) 
 {
  S[i].binRead(is);   
  S[i].setnorm();
 } 
}


void Superellipsoid_D::seta1(double a1) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

void Superellipsoid_D::seta2(double a2) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

void Superellipsoid_D::seta3(double a3) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

void Superellipsoid_D::sete1(double e1) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

void Superellipsoid_D::sete2(double e2) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2); 
}

void Superellipsoid_D::setntheta(double ntheta) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

void Superellipsoid_D::setnphi(double nphi) 
{
  generate(ntheta,nphi,a1,a2,a3,e1,e2);
}

/*-------------------------------------------------------------------------*/
// Hier beginnt die Klasse "Superellipsoid"!!


Superellipsoid::Superellipsoid()
 : Form()
{
 type=SUPERELLIPSOID;
 r0=1.0;
 sf=1.0;
}

Superellipsoid::Superellipsoid (const Superellipsoid &Su):Form(Su)
{
  sf=Su.sf;
  r0=Su.r0;
  sf=1.0;  
  P=Su.P;  
  a1=Su.a1;
  a2=Su.a2;
  a3=Su.a3;
  e1=Su.e1;
  e2=Su.e2;
  initQuad();
  type=SUPERELLIPSOID;
}
Superellipsoid::Superellipsoid (
             const Vector<double> &P,
             double a1, double a2, double a3, double e1, double e2,
             complex<double>  n,
             double r0,
             const Matrix<complex<double> > alpha,
             const Vector<double> &Ex,
             const Vector<double> &Ey,
             const Vector<double> &Ez
             )
: Form (P, n, alpha, Ex, Ey, Ez,SUPERELLIPSOID)
{
 this->a1=a1;
 this->a2=a2;
 this->a3=a3;
 this->e1=e1;
 this->e2=e2;
 this->r0=r0;
 initQuad();
 type=SUPERELLIPSOID;
}

//berechneG berechnet die Lipschitz-Konstante für g(t) in [t1,t2]

double Superellipsoid::berechneG(double &t1, double &t2, Vector<double> p, Vector<double> k)
{
 double bx, by, bz, dgdt;
 double bxs,bx1,bx2,bys,by1,by2,bzs,bz2;
 double A1,A2;
 Vector<double> P1, P2,r;
 P1=p+t1*k;
 P2=p+t2*k;
  r=emax(P1,P2);
 bx=fabs(r[0]/a1);
 by=fabs(r[1]/a2);
 bz=fabs(r[2]/a3);

 if (bx==0.0) { bxs=0.0; bx1=0; bx2=0; }
 else { bxs=pow(bx,2.0/e2); bx1=bxs/bx; bx2=bx1/bx;  }
 
 if (by==0.0) { bys=0.0; by1=0; by2=0; }
 else { bys=pow(by,2.0/e2); by1=bys/by; by2=by1/by;  }

 if (bz==0.0) { bz2=0; }
 else { bzs=pow(bz,2.0/e2); bz2=bzs/(bz*bz);  }
 
 if (bxs+bys==0) {A1=0.0; A2=0.0;}
 else {A1=pow(bxs+bys,e2/e1-2.0); A2=pow(bxs+bys,e2/e1-1.0);}

 dgdt=fabs(2.0/e1*((e2/e1-1.0)*A1*2.0/e2*(k[0]/a1*bx1+k[1]/a2*by1)*(k[0]/a1*bx1+k[1]/a2*by1))
      +A2*(2.0/e1-1.0)*(k[0]*k[0]/(a1*a1)*bx2+k[1]*k[1]/(a2*a2)*by2)+k[2]*k[2]/(a3*a3)*(2.0/e1-1.0)*bz2);
 return dgdt;
}


// berechnegtm berechnet g(tm) in [t1,t2]

double Superellipsoid::berechneg(double &t, Vector<double> p, Vector<double> k)
{
 double  gtm, bx, by, bz,bxs,bys,bx1,by1,bz1,A;
 bx=fabs((p[0]+t*k[0])/a1);
 by=fabs((p[1]+t*k[1])/a2);
 bz=fabs((p[2]+t*k[2])/a3);

 if (bx==0.0) {bxs=0.0; bx1=0.0;}
 else { bxs=pow(bx,2.0/e2); bx1=k[0]/a1*bxs/bx; }

 if (by==0.0) {bys=0.0; by1=0.0;}
 else { bys=pow(by,2.0/e2); by1=k[1]/a2*bys/by; }

 if (bz==0.0) {bz1=0.0;}
 else {bz1=k[2]/a3*pow(bz,2.0/e1-1.0);}
 
 if (bxs+bys==0) A=0;
 else A=pow(bxs+bys,e2/e1-1.0);

 gtm=2.0/e1*A*(bx1+by1)+bz1;
 return gtm;
}

double Superellipsoid::berechneFvont(double t, Vector<double> p, Vector<double> k)
{
 double ft, bx, by, bz;
 bx=fabs((p[0]+t*k[0])/a1);
 by=fabs((p[1]+t*k[1])/a2);
 bz=fabs((p[2]+t*k[2])/a3);
//  if ((bx<0) || (by<0) || (bz<0)) return INF;
 ft=pow(pow(bx,2.0/e2)+pow(by,2.0/e2),e2/e1)+pow(bz,2.0/e1);

 return ft-1;
}

double Superellipsoid::berechneFsvont(double t, Vector<double> p, Vector<double> k)
{
 double fst, bx, by, bz, h1, h2, h3;
 bx=fabs((p[0]+t*k[0])/a1);
 by=fabs((p[1]+t*k[1])/a2);
 bz=fabs((p[2]+t*k[2])/a3);
 h1=pow(pow(bx,2.0/e2)+pow(by,2.0/e2),e2/e1-1.0);
 h2=k[0]/a1*pow(bx,2.0/e2-1.0)+k[1]/a2*pow(by,2.0/e2-1.0);
 h3=k[2]/a3*pow(bz,2.0/e1-1.0);
 fst=2.0/e1*(h1*h2+h3);
 return fst;
}

bool Superellipsoid::naeherung (double &t1, double &t2, Vector<double> p, Vector<double> k, int &RecCounter)
{
 double ft1, ft2, G, gt, t1l, t1r, t2l, t2r, tmi, d, eps, xn, xneu;
 bool found, found1, found2;


if (RecCounter>50) 
{ 
//  cout << "Superellipsoid::naeherung: max. Anzahl Rekursionen überschritten !!!!" << endl; 
  return false; 
}
RecCounter++;
 eps=0.0000001;                 //Abbruch-Parameter
 tmi=(t1+t2)/2.0;
 d=(t2-t1)/2.0;
 G=berechneG(t1,t2,p,k);
 gt=fabs(berechneg(tmi,p,k));
 
 
 
// cout << "d=" << d << "    G=" << G << "   gt=" << gt << "   p=" << p << "    k=" << k << "    t1=" << t1 << "     t2=" << t2 << endl;
if ((G==INF) || (gt==INF)) return false;
 if(gt>G*d)
  {
//   cout << "C" << endl;
   ft1=berechneFvont(t1,p,k);
   ft2=berechneFvont(t2,p,k);
//   cout << "ft1="<< ft1 << " ft2=" << ft2 << endl;
   if (ft1*ft2<0)
     {
//      cout << "Es gibt nur eine Nullstelle" << endl;
      return true;
     }
   else    
    {
//     cout << "Es gibt keine Nullstelle" << endl;
     return false;
    }
  }
 else
  {
   //cout << "else" << endl;
   t1l=t1;
   t1r=tmi;
   t2l=tmi;
   t2r=t2;
   if ((tmi-t1l>=eps) && (t2r-tmi>=eps))
   {
    found1=naeherung(t1l,t1r,p,k,RecCounter);
//    cout << "A" << endl;
    //cout << "found1 fertig" << endl;
        if (found1)
     {
      t1=t1l;
      t2=t1r;
 //     cout << "A1" << endl;
      return true;
     }
    else
     {
      found2=naeherung(t2l,t2r,p,k,RecCounter);
   //   cout << "B" << endl;
    //cout << "found2 fertig" << endl;

     if (found2)
      {
       t1=t2l;
       t2=t2r;
//       cout << "B1" << endl;
       return true;
      }
     }

      return found1 || found2;
    }
    else
     return false;
  }
return found;
}

bool Superellipsoid::next(const Vector<double> &p, const Vector<double> &k, Vector<double> &pout, const int inside)
{
  int RecCounter=0;
  bool n;
  double t1, t2, nst;
  t1=1E-10*r0;
  t2=r0;
/*  cout << "=====================================================" << endl;
  cout << "t1=" << t1 << "     t2=" << t2 << "          p=" << p << "    k=" << k << endl;*/
  n=naeherung(t1,t2,p,k,RecCounter);
  if (n)
   {
    nst=berechneNstNV((t1+t2)/2.0,p,k);
    pout=p+nst*k;
   }
  // cout << "================= ENDE ================================" << endl;
 return n;
}



double Superellipsoid::berechneNstNV (double xn, Vector<double> p, Vector<double> k, double epsilon)
{
   double ft, fst, d, xneu;
   do
   {
    ft=berechneFvont(xn,p,k);
    fst=berechneFsvont(xn,p,k);
    xneu=xn-ft/fst;
    d=fabs(xneu-xn);
    xn=xneu;
   }
   while (d>=epsilon);
   return xn;
}

double Superellipsoid::Volume()
{
 double V, B1, B2;
#ifdef __GNUC__
 B1=tgamma(e1/2.0+1)*tgamma(e1)/tgamma(1.5*e1+1);
 B2=tgamma(e2/2.0)*tgamma(e2/2.0)/tgamma(e2);
#else
 B1=0;
 B2=0;
#endif
 return 2.0*a1*a2*a3*e1*e2*B1*B2;
}

Vector<double> Superellipsoid::norm (const Vector<double> &p)
{
 Vector<double> nv;
 nv[0]=e2/e1*pow(pow(p[0]/a1,2.0/e2)+pow(p[1]/a2,2.0/e2),e2/e1-1)*2.0/(e2*a1)*pow(p[0]/a1,2.0/e2-1);
 nv[1]=e2/e1*pow(pow(p[0]/a1,2.0/e2)+pow(p[1]/a2,2.0/e2),e2/e1-1)*2.0/(e2*a2)*pow(p[1]/a2,2.0/e2-1);
 nv[2]=2.0/(e1*a3)*pow(p[2]/a3,2/e1-1);
 return nv/abs(nv);
}

double Superellipsoid::super(double x, double y, double z)
{
 double erg;
 erg=pow(pow(x/a1,2.0/e2)+pow(y/a2,2.0/e2),e2/e1)-pow(z/a3,2.0/e1);
 return erg;
}

double Superellipsoid::super(Vector<double> r)
{
 return super(r[0],r[1],r[2]);
}

bool Superellipsoid::isInside(const Vector<double> &p)
{
 return super(p)<=0;
}

void Superellipsoid::setParms(double a1,double a2, double a3, double e1, double e2)
{
 this->a1=a1;
 this->a2=a2;
 this->a3=a3;
 this->e1=e1;
 this->e2=e2;
 initQuad();
}

void Superellipsoid::scale( double sf)
{
 a1=a1/this->sf*sf;
 a2=a2/this->sf*sf;
 a3=a3/this->sf*sf;
 this->sf=sf;
}

void Superellipsoid::initQuad()
{
 pul=Vector<double> (-a1,-a2,-a3);
 por=Vector<double> (a1,a2,a3);
 pul=pul+P;
 por=por+P;
}