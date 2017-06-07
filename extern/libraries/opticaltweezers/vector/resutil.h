#include <complex>

#ifndef COMPLEX_I
#define COMPLEX_I
//complex<double> I=complex<double>(0.0,1.0);    // imaginaere Einheit
#define I complex<double>(0.0,1.0)
#endif

// #define INF 1.0/0.0
#define INF -10000
#ifndef RESUTIL_H
#define RESUTIL_H 

#define SENKRECHT  1
#define PARALLEL   0
#define MULTI      1
#define ANZ_STRAHLEN 500     // Anzahl der einfallenden Strahlen

#define IS_NOTHING 0
#define IS_ROW 1
#define IS_COL 2

#include "vector.h"
#include "matrix.h"

#define RAMAN 0
#define FLUORESZENZ 1

#define IN_INC  1 
#define IN_HOST 2
#define IN_INC_AND_HOST 3

typedef struct 
{
 int x,y;
} Point;

typedef struct 
{
 char Name[255];
 double n;
} OptProp;


typedef struct 
{
 /*
    E    : Richtungsvektor des Felds
    EAmp : Amplitude
    k    : Ausbreitungsvektor (normiert)
    phi  : Phase
    pol  : Polarisationsrichtung (senkrecht oder parallel)
    alpha : Winkel von einem Reflexionspunkt zur naechsten  (vom
            Partikelmittelpkt. gerechnet
    m    : Anzahl der bisher durchgefuehrten Strahlumlaeufe
     */


 bool getunnelt;
 int pol;
 Vector<double> k,EAmp;
 complex<double> E;
 complex<double> phi;
 int m; 
 double b; 
} StrahlInfo;

typedef struct 
{
 
 Vector <double> P;
 StrahlInfo      S; 
} StrahlArray;

typedef struct 
{
 /*
    P  : Ort des Einschlusses 
    n  : Brechungsindex  
    a  : Radius (in Einheiten des Partikelradius´)
    alpha : Polarisierbarkeit;
 */
 Vector <double> P;
 complex<double> n;
 double a;
 Matrix<double> alpha;  
} EinschlussInfo;

typedef struct 
{
 int Ebene;
 int Pol;
 int Strahlungsart;
 double angmin,angmax;
 int nang;
 double wave;
 bool isKoherent;
} RRTParmsInfo;



typedef struct 
{
 /*
  nx,ny     : Anzahl der Gitterpunkte in x/y-Richtung
  alpha     : Einfallswinkel (relativ zur x-Achse)
  AnzReflex : Anzahl der max. durch zufuehrenden Reflexionen an der Oberflaeche
              pro Strahl 
  AnzRays   : Anzahl Strahlen
  dx, dy    : Breite und Hoehe einer Gitterzelle 
  dxy       : Diagonale einer Gitterzelle
  r0        : Groesse des Partikels
  l0        : (Vakuum-)Wellenlaenge
  n0        : Brechungsindex des Partikels
  pol       : Polarisation des einfallenden Strahls
  ResRad    : radiale Modenzahl
  ResAzi    : azimutale Modenzahl (phi-Richtung) 
  AnzEin    : Anzahl Einschluesse
  ColMin, ColMax : Minimal-/Maximalwert bei Farbeinteilung (in % des Maximums)
*/
 
 int nx,ny;  
 double alpha,db,dx,dy,dxy; 
 int AnzReflex,AnzRays;
 double bmax,r0,r0end,l0,k0; 
 complex<double> n0;
 int pol; 
 int ResRad, ResAzi;
 int AnzEin;
 double AngleTol, evan;
 double PolAngle;
 int phase;
 bool logscale;
 bool tunneln;
 double ColMax,ColMin; 
 int EinX;
} GlobalParms;

ostream& operator << (ostream& os, GlobalParms parms);
istream& operator >> (istream& is, GlobalParms parms);

double abs2(double x);

void output (int nx, int ny, Vector<complex<double> > **G);
void init_Strahl (GlobalParms Parms, StrahlArray *Strahl);
void sub_Gitter (GlobalParms parms, Vector<complex<double> > **Erg, 
                Vector<complex<double> > **Gitter1,
                Vector<complex<double> > **Gitter2);
void add_Gitter (GlobalParms parms, Vector<complex<double> > **Erg, 
                Vector<complex<double> > **Gitter1,
                Vector<complex<double> > **Gitter2);
void clear (GlobalParms parms, Vector<complex<double> > **Gitter);
void Delete (int n, Vector<complex<double> > **Gitter);
void copy (StrahlInfo &dest, StrahlInfo src);
double minmax (double a, double b);
Vector<double> kart2sph (Vector<double> v);
Vector<double> sph2kart (Vector<double> v);

Point get_grid_point (const GlobalParms &parms,Vector<double> P);
Vector <double> set_grid_point (const GlobalParms& parms, Point P);
double grad (const GlobalParms& parms, Vector<complex<double> > **Gitter, const
      Vector<double> &P);
//complex<double>  asin(const complex<double> &);
complex<double>  acos(const complex<double> &);
Vector<double> next (const GlobalParms& parms, 
                               const Vector<double>& P0,
                               const  Vector<double>& k);

void minmax (double x, double dx, int &min, int &max);
void checkEinschluss (Vector<double>& anf, const Vector<double>& end,
                      StrahlInfo& S, int AnzEin, EinschlussInfo *Ein,
		      Vector<double>& Ps, int &Index);
void checkEinschluss (double r0, Vector<double>& anf, const Vector<double>& end,
                      StrahlInfo& S, int AnzEin, EinschlussInfo *Ein,
		      Vector<double>& Ps, int &Index);
void checkEinschluss (double r0, Vector<double>& anf, const Vector<double>& end,
                      const Vector<double> k, int AnzEin, EinschlussInfo *Ein,
		      Vector<double>& Ps, int &Index);

Vector<double> nextP (Vector<double> P, Vector<double> k, Vector<double> OK,double rK,bool &found); 
void toString (char *S, EinschlussInfo *E, int i);
ostream& operator << (ostream& os, EinschlussInfo E);
      	
bool operator == (EinschlussInfo a, EinschlussInfo b);
ostream &savebinGlobalParms (ostream &os, GlobalParms parms);
istream & loadbinGlobalParms (istream &os, GlobalParms &parms);
GlobalParms readGlobalParms (bool old , ifstream &is);
GlobalParms readGlobalParms (bool old , ifstream *is);
void writeRRTParms (ofstream &os, RRTParmsInfo erg); 
 RRTParmsInfo readRRTParms (bool old, ifstream *is);
RRTParmsInfo readRRTParms (bool old, ifstream &is);
void readRRTParms (ifstream &is,RRTParmsInfo &erg);
#endif
