#ifndef __GITTER_H__
#define __GITTER_H__

#include <complex>
#include "vector.h"
#include "resutil.h"
using namespace std;

class gitter
{

Vector <complex<double> > ****gitterarray;
Vector <complex<double> > ****newarray();
void delarray();

// Testet allgemein, ob eine Gitterzelle ix, iy, iz innerhalb einer 
// Kugel mit Radius r0 und Mittelpunkt x0, y0, z0 liegt
// wird von in_kugel und in_einschluss benoetigt

int ist_innen(int ix, int iy, int iz,
              double x0, double y0, double z0, double r0);
int ist_innen(int ix, int iy, int iz, Vector<double> P, double r0);

double dx2, dy2, dz2;

public:

 double rP;
 int nxmax,nymax,nzmax;
 double xmax, ymax,  zmax, dx, dy, dz;
 Vector<complex<double> > DUMMY;

  gitter();
  gitter(int nx, int ny, int nz, double xmax, double ymax, double zmax);
  gitter(int nx, int ny, int nz, double r);

 ~gitter();

 Vector<complex<double> >& operator() (int  i, int j, int k);
 Vector<complex<double> >& operator() (Vector<int>& P);

                gitter& operator=  (const gitter& g);

// Initialisiert ein Gitter mit anzein Einschluessen, die in ein
// stehen
 void init_gitter(EinschlussInfo *ein, int anzein);

// Erzeugt leeres Gitter
 void init_gitter();

// Ausgabe des Gitters auf Standardausgabe
 void show_gitter();

// Ausgabe des Gitters in Datei fname
 void show_gitter(char* fname);

// Gibt die Koordinaten der Gitterzelle aus, in der der Punkt mit
// den  Koordinaten x0, y0 und z0 liegt
 Vector<int> gitterpunkt(double x0, double y0, double z0);
 Vector<int> gitterpunkt(Vector<double>& P);
// Testet, ob die Gitterzelle ix, iy, iz innerhalb des Partikels liegt
 int  in_kugel(int ix, int iy, int iz);

// Testet, ob die Gitterzelle ix, iy, iz innerhalb des Einschlusses ein
// liegt
 int  in_einschluss(int ix, int iy, int iz, EinschlussInfo ein);

 void set_parms(int nx, int ny, int nz, int xmax, int ymax, int zmax);
};


#endif /*  __GITTER_H__ */
