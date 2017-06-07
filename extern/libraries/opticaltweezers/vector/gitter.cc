#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "gitter.h"
#include <complex>
#include "goodies.h"

using namespace std;



gitter::gitter()
{
  gitterarray = 0;
  dx = 0.0;  dy = 0.0;  dz = 0.0; 
  dx2 = 0.0; dy2 = 0.0; dz2 = 0.0;
  nxmax = 0; nymax = 0; nzmax = 0;
}

gitter::gitter(int nx, int ny, int nz, double x, double y,
               double z)
{
 xmax = x; ymax = y; zmax = z;
 nxmax = nx; nymax = ny; nzmax = nz;

 dx = xmax/nx; dx2 = dx/2.0;
 dy = ymax/ny; dy2 = dy/2.0;
 dz = zmax/nz; dz2 = dz/2.0;

 gitterarray = newarray();
}


gitter::gitter(int nx, int ny, int nz, double r)
{
 xmax = 2.0*r; ymax = xmax; zmax = ymax;
 nxmax = nx; nymax = ny; nzmax = nz;
 rP=r;

 dx = xmax/nx; dx2 = dx/2.0;
 dy = ymax/ny; dy2 = dy/2.0;
 dz = zmax/nz; dz2 = dz/2.0;

 /*cout << "xmax=" << xmax << ", ymax=" << ymax << ", zmax=" << zmax << "\n";
 cout << "nxmax=" << nxmax << ", nymax=" << nymax << ", nzmax=" << nzmax << "\n";
 cout << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << "\n";
 cout << "dx/2=" << dx2 << ", dy/2=" << dy2 << ", dz/2=" << dz2 << "\n";
*/
 gitterarray = newarray();
}

gitter::~gitter()
{
 delarray();
}

// Wird fuer die Konstruktoren der Gitter-Klasse benötigt,
// allokiert gitterarray 

Vector<complex<double> >& gitter::operator()(int i, int j, int k)
{
if (gitterarray[i][j][k]!=0)
 {
   return *gitterarray[i][j][k];
 }
else
   { 
    DUMMY=Vector<complex<double> > (INF,INF,INF);
   return DUMMY;
   }
}

Vector<complex<double> >& gitter::operator()(Vector<int>& P)
{
if (gitterarray[P[0]][P[1]][P[2]]!=0)
 {
   return *gitterarray[P[0]][P[1]][P[2]];
 }
else
   { 
    DUMMY=Vector<complex<double> > (INF,INF,INF);
   return DUMMY;
   }
}
gitter& gitter::operator=(const gitter& g)
{
 if (this != &g)
 {
  xmax = g.xmax; ymax =g.ymax; zmax = g.zmax;
  nxmax = g.nxmax; nymax = g.nymax; nzmax = g.nzmax;
  dx = g.dx; dy  = g.dy; dz = g.dz;
  dx2 = g.dx2; dy2  = g.dy2; dz2 = g.dz2;
  delete gitterarray;
  gitterarray=newarray();
  for (int ix=0; ix<nxmax; ix++)
  {
   for (int iy=0; iy<nymax; iy++)
   {
    for (int iz=0; iz<nzmax; iz++)
    {
     if (g.gitterarray[ix][iy][iz]!=0)
     {
      gitterarray[ix][iy][iz] = new Vector<complex<double> >;
     *gitterarray[ix][iy][iz] = *g.gitterarray[ix][iy][iz];
     }
     else
      gitterarray[ix][iy][iz] = 0;
    }
   }
  }
 }
 return *this;
}

// Wird fuer die Konstruktoren der Gitter-Klasse benötigt,
// allokiert gitterarray

Vector<complex<double> > ****gitter::newarray()
{
Vector<complex<double> > ****hilfarray;

hilfarray = new Vector<complex<double> >***[nxmax];

 for (int ix = 0; ix < nxmax; ix++)
 {
  hilfarray[ix] = new Vector<complex<double> >**[nymax];
  for (int iy = 0; iy < nymax; iy++)
  {
   hilfarray[ix][iy] = new Vector<complex<double> >*[nzmax];
   for (int iz = 0; iz < nzmax; iz++)
   {
    hilfarray[ix][iy][iz] = 0;
   }
  }
 }
 return hilfarray;
}

void gitter::delarray()
{
 for( int ix=0; ix < nxmax; ix++)
 {
  for( int iy=0; iy < nymax; iy++)
  {
   delete [] gitterarray[ix][iy];
  }
  delete [] gitterarray[ix];
 }
 delete [] gitterarray;
}


int gitter::ist_innen(int ix, int iy, int iz,
                      double x0, double y0, double z0, double r0)
// Basisfunktion zur Bestimmung der Lage eines 
// Gitterpunktes ix, iy, iz innerhal einer        
// Kugel vom Radius r_0 mit Mittelpunkt bei
// x_0, y_0, z_0
{
 if (sqr((2*ix+1)*dx2-x0)+sqr((2*iy+1)*dy2-y0)+sqr((2*iz+1)*dz2-z0)<=sqr(r0))
 {
  return 1;
  cout << "innen \n";
 }
 else
  return 0;
} 

int gitter::ist_innen(int ix, int iy, int iz, Vector<double> P, double r0)
// Basisfunktion zur Bestimmung der Lage eines 
// Gitterpunktes ix, iy, iz innerhal einer        
// Kugel vom Radius r_0 mit Mittelpunkt bei P
{
 if (sqr((2*ix+1)*dx2-P[0])+sqr((2*iy+1)*dy2-P[1])+sqr((2*iz+1)*dz2-P[2])<=sqr(r0))
 {
  return 1;
  cout << "innen \n";
 }
 else
  return 0;
}

void gitter::init_gitter(EinschlussInfo *ein, int anzein)
// Initialisierung des Gitters mit anzein Einschlüssen
{
 for(int nein=0;nein<anzein;nein++)
 {
//  double x0 = ein[nein].P[0];
//  double y0 = ein[nein].P[1];
//  double z0 = ein[nein].P[2];
//  double r0 = ein[nein].a;
 
  for(int ix=0;ix<nxmax;ix++)
  {
    for(int iy=0;iy<nymax;iy++)
    {
      for(int iz=0;iz<nzmax;iz++)
      {
       if (gitterarray[ix][iy][iz]==0)
       {
        if (in_einschluss(ix,iy,iz,ein[nein]))
        {
          gitterarray[ix][iy][iz] =  new Vector<complex<double> >;
          *gitterarray[ix][iy][iz] = Vector<complex<double> >(1,1,1);
//          cout << "Einschluss bei ix=" << ix << ", iy= " << iy << ", iz= " << iz << "\n";
//          cout << "test:" <<  *gitterarray[ix][iy][iz] << "\n";
        }
        else
          gitterarray[ix][iy][iz] = 0;
       }
     }
    }
  }
 }
}


int gitter::in_kugel(int ix, int iy, int iz)
// Test, ob ein vorgegebener Punkt ix, iy, iz innerhalb des 
// Partikels liegt
{
 double x0 = xmax/2;
 double y0 = ymax/2;
 double z0 = zmax/2;
 double r0 = x0;
 int hilf = ist_innen(ix, iy, iz, x0, y0, z0, r0);
 return hilf;
}

int gitter::in_einschluss(int ix, int iy, int iz, EinschlussInfo ein)
{
 double r0 = ein.a*rP;
 int hilf = ist_innen(ix, iy, iz, ein.P*rP, r0);
 return hilf;
}

void gitter::show_gitter(char *fname)
{
 ofstream os;
 os.open(fname);
 
 for (int ix=0; ix<nxmax; ix++)
 {
  for (int iy=0; iy<nymax; iy++)
  {
   for (int iz=0; iz<nzmax; iz++)
   {
//   cout << "Gitt:" << gitterarray[ix][iy][iz] << "\n";
   if (gitterarray[ix][iy][iz]!=0)
     {
//      cout << "!!!\n";
      os << (2*ix+1)*dx2 << " " << (2*iy+1)*dy2 << " " << (2*iz+1)*dz2 << "\n";
     }
//   else
//     cout << "kein Einschluss" << "\n";
   }
  }
 }
 os.close();
}

void gitter::show_gitter()
{
 for (int ix=0; ix<nxmax; ix++)
 {
  for (int iy=0; iy<nymax; iy++)
  {
   for (int iz=0; iz<nzmax; iz++)
   {
//   cout << "Gitt:" << gitterarray[ix][iy][iz] << "\n";
   if (gitterarray[ix][iy][iz]!=0)
     {
//      cout << "!!!\n";
     cout << (2*ix+1)*dx2 << " " << (2*iy+1)*dy2 << " " << (2*iz+1)*dz2 << "\n";
     }
//   else
//     cout << "kein Einschluss" << "\n";
   }
  }
 }
}

Vector<int> gitter::gitterpunkt(double x0, double y0, double z0)
{
 Vector<int> punkt;

 punkt[0] = (int)(x0/dx);
 punkt[1] = (int)(y0/dy);
 punkt[2] = (int)(z0/dz);

 return punkt;
}

Vector<int> gitter::gitterpunkt(Vector<double> &P)
{
 Vector<int> punkt;

 punkt[0] = (int)(P[0]/dx);
 punkt[1] = (int)(P[1]/dy);
 punkt[2] = (int)(P[2]/dz);

 return punkt;
}
