/***************************************************************************
                          misc.cpp  -  description
                             -------------------
    begin                : Mit Mär 12 2003
    copyright            : (C) 2003 by Thomas Weigel
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

#include "misc.h"
#include <iostream>
using namespace std;

// extern long int deletedItems, createdItems;

void initInc(Form *E)
{
 switch (E->type)
 {
//  case FUNSURF   : ((Funsurf *)E)->initQuad(); break;
  case ELLIPSOID : ((FormEllipsoid *)E)->initQuad(); break;
  case SUPERELLIPSOID :((Superellipsoid *)E)->initQuad(); break;
  case SUPERELLIPSOID_D :
  case SURFACE   : ((surface *)E)->initQuad(); break;
  case ZYLINDER  : ((Zylinder *)E)->initQuad(); break;
  case KREISKEGEL : ((Kegel *)E)->initQuad(); break;
  case KEGELSTUMPF : ((Kegelstumpf *)E)->initQuad(); break;
  case KEGELSTUMPF_HOHL : ((KegelstumpfHohl *)E)->initQuad(); break;
 }
}

void setR0(Form *E, double r0)
{
 switch (E->type)
 {
  case FUNSURF   : ((Funsurf *)E)->setr0(r0); break;
  case ELLIPSOID : ((FormEllipsoid *)E)->setr0(r0); break;
  case ERYTHROCYTE : ((Erythrocyte *)E)->setr0(r0); break;
  case SUPERELLIPSOID_N : ((Superellipsoid_n *)E)->setr0(r0); break;
  case SUPERELLIPSOID_D :
  case SURFACE   : ((surface *)E)->setr0(r0); break;
  case ZYLINDER  : ((Zylinder *)E)->setr0(r0); break;
  case KREISKEGEL :((Kegel *)E)->setr0(r0); break;
  case KEGELSTUMPF :((Kegelstumpf *)E)->setr0(r0); break;
  case KEGELSTUMPF_HOHL :((KegelstumpfHohl *)E)->setr0(r0); break;
 // case COMPOUND  : ((Compound *)E)->setr0(r0); break;
 }
}

void copyInc (Form *&d, Form *s)
{
// createdItems++;
  switch (s->type)
  { 
    case FUNSURF   : d=new Funsurf (*((Funsurf *)s)); break;
   case ELLIPSOID : d=new FormEllipsoid (*((FormEllipsoid *)s)); break;
   case SUPERELLIPSOID : d=new Superellipsoid (*((Superellipsoid *)s)); break;
   case ERYTHROCYTE : d=new Erythrocyte (*((Erythrocyte *)s)); break;
   case SUPERELLIPSOID_N : d=new Superellipsoid_n (*((Superellipsoid_n *)s)); break;
   case SUPERELLIPSOID_D   : d=new Superellipsoid_D(*((Superellipsoid_D *)s)); break;
   case SURFACE   : d=new surface (*((surface *)s)); break;
   case ZYLINDER  : d=new Zylinder (*((Zylinder *)s)); break;
   case KREISKEGEL  : d=new Kegel (*((Kegel *)s)); break;
   case KEGELSTUMPF  : d=new Kegelstumpf (*((Kegelstumpf *)s)); break;
   case KEGELSTUMPF_HOHL  : d=new KegelstumpfHohl (*((KegelstumpfHohl *)s)); break;
   case COMPOUND  : d=new Compound (*((Compound *)s)); break;
  }
}

void copyFormList (Form **&d, Form **s, int anz)
{
/*
  Kopiert die Einschluss-Liste s in die neue Liste d
*/
 
 if ((anz>0) && (s!=0))
 {
  d=(Form **) malloc (sizeof (Form *) * anz);
 for (int i=0; i<anz; i++)
 {
//  createdItems++;
  switch (s[i]->type)
  {
   case FUNSURF   : d[i]=new Funsurf (*((Funsurf *)s[i])); break;
   case ELLIPSOID : d[i]=new FormEllipsoid (*((FormEllipsoid *)s[i])); break;
   case ERYTHROCYTE : d[i]=new Erythrocyte (*((Erythrocyte *)s[i])); break;
  case SUPERELLIPSOID : d[i]=new Superellipsoid (*((Superellipsoid *)s[i])); break;
  case SUPERELLIPSOID_N : d[i]=new Superellipsoid_n (*((Superellipsoid_n *)s[i])); break;
   case SUPERELLIPSOID_D   : d[i]=new Superellipsoid_D(*((Superellipsoid_D *)s[i])); break;
   case SURFACE   : d[i]=new surface (*((surface *)s[i])); break;
   case ZYLINDER  : d[i]=new Zylinder (*((Zylinder *)s[i])); break;
   case KREISKEGEL  : d[i]=new Kegel (*((Kegel *)s[i])); break;
   case KEGELSTUMPF  : d[i]=new Kegelstumpf (*((Kegelstumpf *)s[i])); break;
   case KEGELSTUMPF_HOHL  : d[i]=new KegelstumpfHohl (*((KegelstumpfHohl *)s[i])); break;
   case COMPOUND  : d[i]=new Compound (*((Compound *)s[i])); break;
   case LINSE     : d[i]=new Linse (*((Linse *)s[i])); break;
  }
//   sprintf (d[i]->Beschreibung,"%s",s[i]->Beschreibung); 
 }
 }
}

void binWriteIncList (ofstream &os, Form **E, int anz)
{
 os.write ((char *) &anz, (char) sizeof (anz));
 if (anz>0)
 for (int i=0; i<anz; i++)
  binWriteInc(os,E[i]);
}

void binWriteInc (ofstream &os, Form *E)
{
 os.write ((char *) &E->type, (char) sizeof (E->type));
 switch (E->type)
 {
  case FUNSURF   : ((Funsurf *)E)->binWrite(os); break;
  case ELLIPSOID : ((FormEllipsoid *)E)->binWrite(os); break;
  case SUPERELLIPSOID_D : ((Superellipsoid_D *)E)->binWrite(os); break;
  case SUPERELLIPSOID_N : ((Superellipsoid_n *)E)->binWrite(os); break;
  case ERYTHROCYTE : ((Erythrocyte *)E)->binWrite(os); break;
  case SURFACE   : ((surface *)E)->binWrite(os); break;
  case ZYLINDER  : ((Zylinder *)E)->binWrite(os); break;
  case KREISKEGEL  : ((Kegel *)E)->binWrite(os); break;
  case KEGELSTUMPF  : ((Kegelstumpf *)E)->binWrite(os); break;
  case KEGELSTUMPF_HOHL  : ((KegelstumpfHohl *)E)->binWrite(os); break;
  case COMPOUND  : ((Compound *)E)->binWrite(os); break;
 }
}

void binReadIncList (ifstream &is, Form **&E, int anz)
{
 int Anz;
 is.read ((char *) &Anz, (char) sizeof (Anz));
 cout << "ANZ=" << Anz << endl;
 if (anz>0)
 {
 E=(Form **) malloc (sizeof (Form *) * anz);
 for (int i=0; i<anz; i++)
 {
  binReadInc(is,E[i],true);
 }
 }
 else E=0; 
}

void binReadInc (ifstream &is, Form *&E, bool isNew)
{
 int type;
 is.read ((char *) &type, (char) sizeof (type)); 
// createdItems++;
 switch (type)
 {
  case FUNSURF   : if (isNew) E=new Funsurf; ((Funsurf *)E)->binRead(is); break;
  case ELLIPSOID : if (isNew) E=new FormEllipsoid; ((FormEllipsoid *)E)->binRead(is); break;
  case ERYTHROCYTE : if (isNew) E=new Erythrocyte; ((Erythrocyte *)E)->binRead(is); break;
  case SUPERELLIPSOID_D   : if (isNew) E=new Superellipsoid_D; ((Superellipsoid_D *)E)->binRead(is); break;
  case SUPERELLIPSOID_N   : if (isNew) E=new Superellipsoid_n; ((Superellipsoid_n *)E)->binRead(is); break;
  case SURFACE   : if (isNew) E=new surface; ((surface *)E)->binRead(is); break;
  case ZYLINDER  : if (isNew) E=new Zylinder; ((Zylinder *)E)->binRead(is); break;
  case KREISKEGEL  : if (isNew) E=new Kegel; ((Kegel *)E)->binRead(is); break;
  case KEGELSTUMPF  : ((Kegelstumpf *)E)->binRead(is); break;
  case KEGELSTUMPF_HOHL  : ((KegelstumpfHohl *)E)->binRead(is); break;
  case COMPOUND : if (isNew) E=new Compound; ((Compound *)E)->binRead(is); break;
  }
}

void deleteInc (Form *E)
{
// deletedItems++;
/*
 switch (E->type)
 {
  case FUNSURF : delete ((Funsurf *) E); break;
  case ELLIPSOID : delete ((Ellipsoid *) E); break;
  case SUPERELLIPSOID_D : delete ((Superellipsoid_D *)E); break;
  case SUPERELLIPSOID : delete ((Superellipsoid *)E); break;
  case SURFACE : delete ((surface *)E); break;
  case ZYLINDER : delete ((Zylinder *)E); break;
  case COMPOUND : delete ((Compound *)E); break;
 }*/
 delete E;
 E=0;
}
double NA2w0(double lambda, double NA, complex<double> n)
{
	double nr=real(n);
	double w0=lambda/M_PI*sqrt(nr*nr-NA*NA)/NA; 
	return w0;
}

Vector<double> force (Vector<double> norm, Ray Se, Ray Sr, Ray St, double df)
// hier wird mu_x als 1 angenommen !
// norm zeigt entgegen der Einfallsrichtung des Strahls !!!!
{
   Vector<double> F,fe,fr,ft,n;
   // double dF=Se.flaeche();
   double dF=df;
   double ne=real(Se.n);
   double nr=real(Sr.n);
   double nt=real(St.n);
  // double h=Se.k0/Z0; 
   double h=1.0/c_light;   
   if (!Se.imEinschluss) fe=ne*dF*abs2(Se.E[4])*Se.k[4]*fabs(norm*Se.k[4]);
   if (!Sr.imEinschluss) fr=nr*dF*abs2(Sr.E[4])*Sr.k[4]*fabs(norm*Sr.k[4]);
   if (!St.imEinschluss) ft=nt*dF*abs2(St.E[4])*St.k[4]*fabs(norm*St.k[4]); 

   F=h*(fe-fr-ft);
 
 //  cout << Se.P[4] << "    " << F << endl;
   return F;
}

/*Vector<double> force (Vector<double> norm, Strahl Se, Strahl Sr, Strahl St, double df)
{
	double dA=Se.flaeche();
	double phi_r,phi_t;
    Vector<double> F;
	phi_r=0.5*c_light*real(Sr.n)*eps0*dA*abs2(Sr.E[4]);
	phi_t=0.5*c_light*real(St.n)*eps0*dA*abs2(St.E[4]);
	F=real(Se.n)/c_light*((Sr.k[4]-Se.k[4])*phi_r+(St.k[4]-Se.k[4])*phi_t);
	return F;
}*/

double gaussw(double z, double wvl, double w0)
{
 double z0=M_PI*w0*w0/wvl;
 double w=w0*sqrt(1+z*z/(z0*z0));
 return w;
}

complex<double> gaussphase (Vector<double> P, Vector<double> F, Vector<double> k, double w0, double k0)
// Achtung k wird als normiert angenommen !
{
	Vector<double> h=F-P;
	double z=k*h;
	double r=abs(h-z*k);
    double z0=k0*w0/2.0;
	double rz=z/z0;
	double zeta=atan(rz);
	double R=z*(1+rz*rz);
	return exp(-I*k0*r*r/2.0/R)*exp(I*(zeta-k0*z));   
}

float readLE_float32(istream &is)
{
	char *d;
	char h;
	float f;
	d=(char *)&f;
	is.read (d,4);
	/*h=d[3];
	d[3]=d[0];
	d[0]=h;
	h=d[1];
	d[1]=d[2];
	d[2]=h;*/
	return f;
}

int readLE_int32(istream &is)
{
	unsigned char d[4];
	is.read((char *)d,4);
	int i=d[0]+d[1]*256+d[2]*65536+d[3]*16777216;
	return i;
}