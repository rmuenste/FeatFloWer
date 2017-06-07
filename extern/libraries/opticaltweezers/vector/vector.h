#ifndef VEKTOR_H
#define VEKTOR_H
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <fstream>
#define _USE_MATH_DEFINES
#include<math.h> 

using namespace std; 
/*
#ifndef MIN_XY
#define MIN_XY
#define min(x,y) x<y ? x : y
#endif 

#ifndef MAX_XY
#define MAX_XY
#define max(x,y) x>y ? x : y
#endif 

*/
 
double arctan (double y, double x);
complex<double>  tan(complex<double>  z);

inline complex<double> operator + (complex<double> x, int y) { return x+(double)y; }
inline complex<double> operator + (int y, complex<double> x) { return (double)y+x; }
inline complex<double> operator - (complex<double> x, int y) { return x-(double)y; }
inline complex<double> operator - (int y, complex<double> x) { return (double)y-x; }

inline complex<double> operator * (complex<double> x, int y) { return x*(double)y; }
inline complex<double> operator * (int y, complex<double> x) { return (double)y*x; }
inline complex<double> operator / (complex<double> x, int y) { return x/(double)y; }
inline complex<double> operator / (int y, complex<double> x) { return (double)y/x; }

/** Template-Klasse für dreidimensionale Vektoren
*/
template <class T>
class Vector            // Das ist die eigentliche Vektor-Klasse
{
 public :
 //! Binäres Schreiben eines Vektors (Komponenten werden nacheinander binär gespeichert)
 void binWrite (ofstream &os)
 {
  for (int i=0; i<3; i++)
  os.write ((char *) &data[i], (char) sizeof (data[i]));
 }

//! Binäres Lesen eines Vektors 
 void binRead (ifstream &is)
 {
  for (int i=0; i<3; i++)
  {
   is.read ((char *) &data[i], (char) sizeof (data[i]));   
  }  
 }
 /** Klammer-Operatoren: 
   Ohne Argument werden die Vektoren auf (0,0,0) initalisiert */
 Vector () 
  { 
   data[0]=0;
   data[1]=0;
   data[2]=0;
  }

 /** Klammer-Operator mit den einzelnen Komponenten als Argument */
 Vector (T x, T y, T z) 
  { 
   data[0]=x; 
   data[1]=y;
   data[2]=z;
  }
 
/** Copy-Konstruktor */ 
inline Vector (const Vector& r)
 { 
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
 }
  
/** Addition Vektor mit Vektor */
inline  friend Vector operator + (const Vector& r1, const Vector& r2)
 {
 Vector<T> h; 
  h.data[0]=r1.data[0] + r2.data[0];
  h.data[1]=r1.data[1] + r2.data[1];
  h.data[2]=r1.data[2] + r2.data[2];
 return h;
 }

/** Vorzeichen-Minus */
inline Vector operator - ()
 {
  Vector Erg;
  for (int i=0; i<3; i++)
  Erg.data[i]=-data[i];
  return Erg;
 } 
 
/** Subtraktion Vektor mit Vektor */
inline friend Vector operator - (const Vector& r1, const Vector& r2)
 {
  Vector<T> h;
  h.data[0]=r1.data[0] - r2.data[0];
  h.data[1]=r1.data[1] - r2.data[1];
  h.data[2]=r1.data[2] - r2.data[2];
  return h; 
 }
 

/** Kreuzprodukt */ 
inline friend Vector operator % (const Vector &r1,const Vector& r2)
 // Kreuzprodukt
 {
 Vector<T> h;
 h.data[0]=r1.data[1]*r2.data[2] - r1.data[2]*r2.data[1];
 h.data[1]=r1.data[2]*r2.data[0] - r1.data[0]*r2.data[2];
 h.data[2]=r1.data[0]*r2.data[1] - r1.data[1]*r2.data[0];
 return h;
 }

/** Quadratwurzel aus einem Vektor: Ergebnis: sqrt((x,y,z)) -> (sqrt(x),sqrt(y),sqrt(z)) */
inline friend Vector sqrt (const Vector& r)
 {
  return Vector<T> (sqrt(r.data[0]),sqrt(r.data[1]),sqrt(r.data[2]));
 }

/** Exponentialfunktion Ergebnis: exp((x,y,z)) -> (exp(x),exp(y),exp(z)) */
inline  friend Vector exp(const Vector& r)
 {
  return Vector<T> (exp(r.data[0]),exp(r.data[1]),exp(r.data[2]));
 }

/** Komponentenweise Multiplikation zweier Vektoren */ 
inline friend Vector emult (const Vector& r1, const Vector& r2)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)  
    Erg.data[i]=r1.data[i]*r2.data[i];
  return Erg;  
 } 

/** Komponentenweise Division zweier Vektoren */
inline friend Vector ediv (const Vector& r1, const Vector& r2)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)  
    Erg.data[i]=r1.data[i]/r2.data[i];
  return Erg;  
 } 


/* Vector operator * (T x)
 {
  int i;
  Vector<T> h;
  for (i=0; i<3; i++)
   h[i]=data[i] * x;
  return h;
 }*/
 
/** Zuweisungsoperator */ 
 Vector& operator = (const Vector& r)
 {
  if (this == &r) return *this;
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
  return *this;
 }

/** Zuweisungsoperator */ 
 Vector& operator = (T r[3])
 {
  for (int i=0; i<3; i++)
  data[i]=r[i];
  return *this;
 }

/** Vergleichsoperatoren: Zwei Vektoren sind dann gleich, wenn ihre Komponenten übereinstimmen */
  bool operator == (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]==r.data[0]) && (data[1]==r.data[1]) && (data[2]==r.data[2]));
  return Erg;
 } 
 
/** Ungleichoperator: Zwei Vektoren sind ungleich, wenn eine Komponente nicht übereinstimmt */
 bool operator != (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]!=r.data[0]) || (data[1]!=r.data[1]) || (data[2]!=r.data[2]));
  return Erg;
 }
 
/** Addition */
 Vector& operator += (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]+r.data[i];
  return *this;
 }
 
/** Subtraktion */
 Vector& operator -= (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]-r.data[i];
  return *this;
 }
 
/** Multiplikation mit Skalar */
 Vector& operator *= (T r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]*r;
  return *this;
 }
 
/** Division mit Skalar */
 Vector& operator /= (T r) 
 {
   //if (abs(r)==0.0) errmsg ("Operator /: Division durch Null");
   for (int i=0; i<3; i++)
   data[i]=data[i]/r;
   return *this;
 }
 
inline friend bool operator == (const Vector& A, const Vector& B)
 {
  return ((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }
 
inline friend bool operator != (const Vector& A, const Vector& B)
 {
 return !((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }


/** Skalarprodukt Vektor * Vektor */ 
inline friend T operator * (const Vector& a, const Vector& b)
 {
  T h;
  h=0;
  for (int i=0; i<3; i++)
   h+=a[i]*b[i];
  return h;
 }

/** Multiplikation Skalar mit Vektor */
 friend Vector operator * (T x, const Vector& r)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }

/** Multiplikation Vektor mit Skalar */
 friend Vector operator * ( const Vector& r, T x)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }
 

/** Division durch Skalar */
 friend Vector operator / (const Vector& r,T x)
 {
  Vector<T> Erg;
  //if (x==0) errmsg ("Operator /: Division durch Null");
  for (int i=0; i<3; i++)
   Erg.data[i]=r.data[i] / x;
  return Erg; 
 }
 

/** Zugriff auf einzelne Komponente (wie bei normalem Array) */
 T& operator[] (int i) {return data[i]; };

 const T& operator[ ] (int i) const  { return data[i]; };
 T data[3]; 
 //int Dim;
};

template <class T> void out_vector (int Anz,Vector<T> *r)
{
 for (int i=0; i<Anz; i++)
 cout << r[i] << endl; 
}


template <class T> Vector<T> cart2sphere(const Vector<T>& x)

/**
  transformiert kartesische Koordinaten in Kugelkoordinaten

  Eingang: x  (kart. Koord.)

  Ausgang:   (Kugel-Koord.)
*/

{
 Vector<T> y,a;

  a[0] = x[0];
  a[1] = x[1];
  a[2] = x[2];
  y[1] = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  if (a[2] == 0.0)
    if (y[0] == 0.0)
      y[1] = 0;
    else
      y[1] = M_PI / 2.0;
    else
      {
        y[1] = atan(sqrt(a[0] * a[0] + a[1] * a[1]) / a[2]);
        if (a[2] < 0.0)
          y[1] = M_PI + y[1];
      }
  if (a[0] == 0.0)
    if (a[1] == 0.0)
      y[2] = 0.0;
    else
      if (a[1] > 0.0)
        y[2] = M_PI / 2.0;
      else
        y[2] = 3.0 * M_PI / 2.0;
  else
    {
      y[2] = atan(a[1] / a[0]);
      if (a[0] < 0.0)
        y[2] = M_PI + y[2];
    }
  return y;
}




template <class T> Vector<T> sphere2cart (const Vector<T>& P)
{
 /** Wandelt den sphaerischen Vektor P in einen kartesischen Vektor um
     gemaess : (r,theta,phi) -> (x,y,z) */

 Vector<T> Erg;
 Erg[0]=P[0]*sin(P[1])*cos(P[2]);
 Erg[1]=P[0]*sin(P[1])*sin(P[2]);
 Erg[2]=P[0]*cos(P[1]); 
 return Erg; 
}


istream& operator >> (istream &is, Vector<complex<double> > &r);
istream& operator >> (istream &is, Vector<double> &r);
istream& operator >> (istream &is, Vector<int> &r);

ostream& operator << (ostream &is, const Vector<complex<double> > &r);
ostream& operator << (ostream &is, const Vector<double> &r);
ostream& operator << (ostream &is, const Vector<int> &r);

void errmsg(char *S); //! Fehlermeldung ausgeben
double abs (const Vector<double> &r); //! Betrag eines double-Vektors
double abs (const Vector<complex<double> > &r); //! Betrag eines komplexen Vektors
double abs (const Vector<int> &r); //! Betrag eines integer-Vektors
double abs2 (const Vector<double> &r); //! Betragsquadrat eines double-Vektors
double abs2 (const Vector<complex<double> > &r);//! Betragsquadrat eines komplexen Vektors
double  abs2 (complex<double>  x); //! Betragsquadrat einer komplexen Zahl
double abs2 (double x); //! Betragsquadrat einer double Zahl

Vector<double> operator - (const Vector<double>&r1, const Vector<int> &r2);
Vector<double> operator - (const Vector<int>&r1, const Vector<double> &r2);
Vector<complex<double> > operator - (const Vector<double> &r1,const Vector<complex<double> > &r2);
Vector<complex<double> > operator - (const Vector<complex<double> > &r1,const Vector<double> &r2);
Vector<complex<double> > operator - (const Vector<complex<double> >& r1, const Vector<int> &r2);
Vector<complex<double> > operator - (const Vector<int>& r1, const Vector<complex<double> > &r2);

Vector<complex<double> > operator + (const Vector<double> &r1,const Vector<complex<double> > &r2);
Vector<complex<double> > operator + (const Vector<complex<double> > &r1,const Vector<double> &r2);
Vector<double> operator + (const Vector<int>& r1, const Vector<double> &r2);
Vector<double> operator + (const Vector<double>& r1, const Vector<int> &r2);
Vector<complex<double> > operator + (const Vector<complex<double> >& r1, const Vector<int> &r2);
Vector<complex<double> > operator + (const Vector<int>& r1, const Vector<complex<double> > &r2);

Vector<complex<double> > conj (const Vector<complex<double> > &r);

Vector<double> operator * (int x, const Vector<double> &r);
Vector<double> operator * (const Vector<double> &r,int x);
Vector<double> operator * (double x,const Vector<int>& r);
Vector<double> operator * (const Vector<int>& r, double x);

/*Vector<complex<double> > operator * (double x,const Vector<complex<double> >& r);
Vector<complex<double> > operator * (const Vector<complex<double> >& r,double x);
Vector<complex<double> > operator * (int x,const Vector<complex<double> >& r);
Vector<complex<double> > operator * (const Vector<complex<double> >& r,int x);
*/

Vector<complex<double> > operator * (complex<double>  x,const Vector<double>& r);
Vector<complex<double> > operator * (const Vector<double>& r, complex<double>  x);
Vector<complex<double> > operator * (complex<double>  x,const Vector<int>& r);
Vector<complex<double> > operator * (const Vector<int>& r, complex<double>  x);

Vector<complex<double> > operator / (const Vector<double>& r, complex<double>  x);
Vector<double> operator / (const Vector<double>& r, int x);
Vector<complex<double> > operator / (const Vector<complex<double> >& r, double x);
Vector<complex<double> > operator / (const Vector<complex<double> >& r, int x);
Vector<double> operator / (const Vector<int>& r, double x);
Vector<complex<double> > operator / (const Vector<int>& r, complex<double>  x);


complex<double>  operator * (const Vector<double>& a, const Vector<complex<double> > &b);
complex<double>  operator * (const Vector<complex<double> >& a, const Vector<double> &b);
double operator * (const Vector<double>& a, const Vector<int> &b);
double operator * (const Vector<int>& a, const Vector<double> &b);
complex<double>  operator * (const Vector<complex<double> >& a, const Vector<int> &b);
complex<double>  operator * (const Vector<int>& a, const Vector<complex<double> > &b);

Vector<complex<double> > operator % (const Vector<complex<double> > &a, 
                                   const Vector<double> &b);
Vector<complex<double> > operator % (const Vector<double> &a, 
                                   const Vector<complex<double> > &b);
Vector<complex<double> > operator % (const Vector<complex<double> > &a, const Vector<int> &b);
Vector<complex<double> > operator % (const Vector<int> &a, const Vector<complex<double> > &b);
Vector<complex<double> > operator % (const Vector<double> &a, const Vector<int> &b);
Vector<complex<double> > operator % (const Vector<int> &a, const Vector<double> &b);


double *conv2double (int Anz, Vector<complex<double> > *r);

Vector<double> emult (const Vector<double> &r1, const Vector<int> &r2);
Vector<double> emult (const Vector<int> &r1, const Vector<double> &r2);
Vector<double> floor (const Vector<double> &r);
Vector<int> ifloor (const Vector<double> &r);
inline Vector <double> ceil (const Vector<double>& r)
{
 return Vector<double> ((int) ceil (r[0]),(int) ceil (r[1]), (int) ceil (r[2]));
}

inline Vector <int> iceil (const Vector<double>& r)
{
 return Vector<int> ((int) ceil (r[0]),(int) ceil (r[1]), (int) ceil (r[2]));
}

inline double sign(double x)
{
 double y=(x>0) ? 1.0 : 0.0 + (x<0) ? -1.0 : 0.0;
 return y;
}
Vector<complex<double> > *conv2vector (int Anz, double *r);
Vector<complex<double> > grad2d (double (*f(Vector<complex<double> >)),Vector<complex<double> > x, double dx);
Vector<complex<double> > convd2c (const Vector<double> &r); 
Vector<double> real(Vector <complex<double> > r);
Vector<double> imag(Vector <complex<double> > r);
double sDreieck (Vector<double>,Vector<double>,Vector<double>);
double sViereck (Vector<double>,Vector<double>,Vector<double>,Vector<double>); 
void dreh (Vector<double> &r, double phi);
Vector<double> dreh (const Vector<double> &r, double dtheta, double dphi);
Vector<complex<double> > vdc (const Vector<double> &r);
Vector<complex<double> > makeReal (const Vector<complex<double> > &r);
Vector<double> cart2sphere(double x, double y, double z);
Vector<double> sphere2cart (double r, double theta, double phi);
void getKSystem (const Vector<double> &n, const Vector<double> &k,
                Vector<double> &e1, Vector<double> &e2, Vector<double> &e3);

Vector<double> arg (Vector<complex<double> > &r);
complex<double>  asin (complex<double>  z);
complex<double>  tan  (complex<double>  z);
complex<double>  ihoch (int l);
void sphereunitv (Vector<double> &P, Vector<double> &er, Vector<double> &etheta, Vector<double> &ephi); 
Vector <double> emax (Vector<double> &P1, Vector<double> &P2);
/*const Vector<double> ex=Vector<double> (1.0,0.0,0.0);
const Vector<double> ey=Vector<double> (0.0,1.0,0.0);
const Vector<double> ez=Vector<double> (0.0,0.0,1.0);
const Vector<double> zero=Vector<double> (0.0,0.0,0.0);
const Vector<double>one=Vector<double>(1.0,1.0,1.0);*/
#define ex Vector<double> (1.0,0.0,0.0)
#define ey Vector<double> (0.0,1.0,0.0)
#define ez Vector<double> (0.0,0.0,1.0)
#define zero Vector<double> (0.0,0.0,0.0)
#define one Vector<double>  (1.0,1.0,1.0)
#define czero Vector<complex<double> > (0.0,0.0,0.0)
#define cone Vector<complex<double> >  (1.0,1.0,1.0)

Vector<double> er(double theta, double phi);
Vector<double> etheta(double theta, double phi);
Vector<double> ephi(double theta, double phi);

Vector<complex<double> > operator * (complex<double>  x,const Vector<double>& r);
Vector<complex<double> > operator * (const Vector<double>& r, complex<double>  x);
double *conv2double (int Anz, Vector<complex<double> > *r);
Vector<complex<double> > *conv2vector (int Anz, double *r);
Vector<complex<double> > grad2d (double (*f(Vector<complex<double> >)),Vector<complex<double> > x, double dx);
Vector<complex<double> > convd2c (const Vector<complex<double> > &r); 
void dreh (Vector<double> &r, double phi);
complex<double>  asin (complex<double>  z);
complex<double>  tan  (complex<double>  z);
Vector<complex<double> > conj (const Vector<complex<double> > &r);
Vector<complex<double> > norm (const Vector<complex<double> > &r);  // Normiert Vektor r
Vector<double> norm (const Vector<double> &r);  // Normiert Vektor r (r=r/|r|)

typedef struct
{
 double a,b,c;             // Seitenlaengen
 double alpha,beta,gamma;  // Winkel
} STriangle;

/* template <class T> Vector<T> sphere2cart (const double& r, const double& theta, const double& phi)
{
 Vector<T> Erg;
 Erg[0]=r * cos(phi) * sin(theta);
 Erg[1]=r * sin(phi) * sin(theta);
 Erg[2]=r * cos(theta);
 return Erg;
}*/


char *toString (char *,Vector<double>);
char* toString (char *S, Vector<complex<double> > P);
 Vector<double> nanV (const char *tagp);   /// gibt einen NaN-Vektor zurück


 bool isnan (Vector<double> v); /// gibt true zurück, wenn eine der Komponenten NaN ist

 #endif

#ifndef COMPLEX_I
#define COMPLEX_I
const complex<double>  I=complex<double>  (0.0,1.0);
#endif


