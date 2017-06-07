#ifndef VEKTOR_H
#define VEKTOR_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <fstream>

#ifndef double_complex
#define double_complex complex<double>
#endif

using namespace std;

/*#ifndef min(x,y)
#define min(x,y) x<y ? x : y
#endif 

#ifndef max(x,y)
#define max(x,y) x>y ? x : y
#endif 
*/


double arctan (double y, double x);
double_complex tan(double_complex z);

inline complex<double> operator + (complex<double> x, int y) { return x+(double)y; }
inline complex<double> operator + (int y, complex<double> x) { return (double)y+x; }
inline complex<double> operator - (complex<double> x, int y) { return x-(double)y; }
inline complex<double> operator - (int y, complex<double> x) { return (double)y-x; }

inline complex<double> operator * (complex<double> x, int y) { return x*(double)y; }
inline complex<double> operator * (int y, complex<double> x) { return (double)y*x; }
inline complex<double> operator / (complex<double> x, int y) { return x/(double)y; }
inline complex<double> operator / (int y, complex<double> x) { return (double)y/x; }


template <class T>
class Vector            // Das ist die eigentliche Vektor-Klasse
{
 public :

 void binWrite (ofstream &os)
 {
  for (int i=0; i<3; i++)
  os.write ((char *) &data[i], (char) sizeof (data[i]));
 }

 void binRead (ifstream &is)
 {
  for (int i=0; i<3; i++)
  is.read ((char *) &data[i], (char) sizeof (data[i]));
 }
  
 Vector () 
  { 
   data[0]=0;
   data[1]=0;
   data[2]=0;
  }

 Vector (T x, T y, T z) 
  { 
   data[0]=x;
   data[1]=y;
   data[2]=z;
  }
  
inline Vector (const Vector& r)
 { 
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
 }
  
inline  friend Vector operator + (const Vector& r1, const Vector& r2)
 {
 Vector<T> h;
  h.data[0]=r1.data[0] + r2.data[0];
  h.data[1]=r1.data[1] + r2.data[1];
  h.data[2]=r1.data[2] + r2.data[2];
 return h;
 }

 Vector operator - ()
 {
  Vector Erg;
  for (int i=0; i<3; i++)
  Erg.data[i]=-data[i];
  return Erg;
 } 
 
inline friend Vector operator - (const Vector& r1, const Vector& r2)
 {
  Vector<T> h;
  h.data[0]=r1.data[0] - r2.data[0];
  h.data[1]=r1.data[1] - r2.data[1];
  h.data[2]=r1.data[2] - r2.data[2];
  return h; 
 }
 
 
 friend Vector operator % (const Vector &r1,const Vector& r2)
 // Kreuzprodukt
 {
 Vector<T> h;
 h.data[0]=r1.data[1]*r2.data[2] - r1.data[2]*r2.data[1];
 h.data[1]=r1.data[2]*r2.data[0] - r1.data[0]*r2.data[2];
 h.data[2]=r1.data[0]*r2.data[1] - r1.data[1]*r2.data[0];
 return h;
 }

inline friend Vector sqrt (const Vector& r)
 {
  return Vector<T> (sqrt(r.data[0]),sqrt(r.data[1]),sqrt(r.data[2]));
 }

inline  friend Vector exp(const Vector& r)
 {
  return Vector<T> (exp(r.data[0]),exp(r.data[1]),exp(r.data[2]));
 }

 
inline friend Vector emult (const Vector& r1, const Vector& r2)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++)  
    Erg.data[i]=r1.data[i]*r2.data[i];
  return Erg;  
 } 

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
  
 Vector& operator = (const Vector& r)
 {
  if (this == &r) return *this;
  for (int i=0; i<3; i++)
  data[i]=r.data[i];
  return *this;
 }
 
  bool operator == (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]==r.data[0]) && (data[1]==r.data[1]) && (data[2]==r.data[2]));
  return Erg;
 } 
 
 bool operator != (const Vector& r)
 {
  bool Erg;
  Erg=((data[0]!=r.data[0]) || (data[1]!=r.data[1]) || (data[2]!=r.data[2]));
  return Erg;
 }
 
 Vector& operator += (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]+r.data[i];
  return *this;
 }
 
 Vector& operator -= (const Vector& r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]-r.data[i];
  return *this;
 }
 
 Vector& operator *= (T r)
 {
  for (int i=0; i<3; i++)
  data[i]=data[i]*r;
  return *this;
 }
 
 Vector& operator /= (T r) 
 {
   //if (abs(r)==0.0) errmsg ("Operator /: Division durch Null");
   for (int i=0; i<3; i++)
   data[i]=data[i]/r;
   return *this;
 }
 
 friend bool operator == (const Vector& A, const Vector& B)
 {
  return ((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }
 
 friend bool operator != (const Vector& A, const Vector& B)
 {
 return !((A[0]==B[0])&&(A[1]==B[1])&&(A[2]==B[2])); 
 }
 
 friend T operator * (const Vector& a, const Vector& b)
 {
  T h;
  h=0;
  for (int i=0; i<3; i++)
   h+=a[i]*b[i];
  return h;
 }

 friend Vector operator * (T x, const Vector& r)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }

 friend Vector operator * ( const Vector& r, T x)
 {
  Vector<T> h;
  for (int i=0; i<3; i++)
  h.data[i]=r.data[i]*x;
  return h; 
 }
 
 friend Vector operator / (const Vector& r,T x)
 {
  Vector<T> Erg;
  //if (x==0) errmsg ("Operator /: Division durch Null");
  for (int i=0; i<3; i++)
   Erg.data[i]=r.data[i] / x;
  return Erg; 
 }
 

 T& operator[] (int i) {return data[i]; };

 const T& operator[ ] (int i) const  { return data[i]; };
 //int dim() { return 3; }

 T data[3]; 
 //int Dim;
};


template <class T> void out_vector (int Anz,Vector<T> *r)
{
 for (int i=0; i<Anz; i++)
 cout << r[i] << endl; 
}

template <class T> Vector<T> cart2sphere (const Vector<T>& P)
{
 //  Wandelt den kartesischen Vektor P in einen sphaerischen Vektor um 
 //  gemaess (x,y,z) -> (r,theta,phi)

 Vector<T> Erg;
 Erg[0]=abs(P);
 if (P[2]!=0.0) Erg[1]=atan2(sqrt(P[0]*P[0]+P[1]*P[1]),P[2]); else Erg[1]=0.0;
 if (P[0]!=0.0) Erg[2]=atan2(P[1],P[0]); else Erg[2]=0.0;
 if (Erg[2]<0.0) Erg[2]=Erg[2]+2.0*M_PI;
 if (Erg[1]<0.0) Erg[1]=Erg[1]+2.0*M_PI;
 return Erg;  
}

template <class T> Vector<T> sphere2cart (const Vector<T>& P)
{
 // Wandelt den sphaerischen Vektor P in einen kartesischen Vektor um
 // gemaess : (r,theta,phi) -> (x,y,z)

 Vector<T> Erg;
 Erg[0]=P[0]*sin(P[1])*cos(P[2]);
 Erg[1]=P[0]*sin(P[1])*sin(P[2]);
 Erg[2]=P[0]*cos(P[1]); 
 return Erg; 
}


istream& operator >> (istream &is, Vector<double_complex> &r);
istream& operator >> (istream &is, Vector<double> &r);
istream& operator >> (istream &is, Vector<int> &r);

ostream& operator << (ostream &is, const Vector<double_complex> &r);
ostream& operator << (ostream &is, const Vector<double> &r);
ostream& operator << (ostream &is, const Vector<int> &r);

void errmsg(char *S);
double abs (const Vector<double> &r);
double abs (const Vector<double_complex> &r);
double abs (const Vector<int> &r);
double abs2 (const Vector<double> &r);
double abs2 (const Vector<double_complex> &r);
double_complex abs2 (double_complex x);
double abs2 (double x);

Vector<double> operator - (const Vector<double>&r1, const Vector<int> &r2);
Vector<double> operator - (const Vector<int>&r1, const Vector<double> &r2);
Vector<double_complex> operator - (const Vector<double> &r1,const Vector<double_complex> &r2);
Vector<double_complex> operator - (const Vector<double_complex> &r1,const Vector<double> &r2);
Vector<double_complex> operator - (const Vector<double_complex>& r1, const Vector<int> &r2);
Vector<double_complex> operator - (const Vector<int>& r1, const Vector<double_complex> &r2);

Vector<double_complex> operator + (const Vector<double> &r1,const Vector<double_complex> &r2);
Vector<double_complex> operator + (const Vector<double_complex> &r1,const Vector<double> &r2);
Vector<double> operator + (const Vector<int>& r1, const Vector<double> &r2);
Vector<double> operator + (const Vector<double>& r1, const Vector<int> &r2);
Vector<double_complex> operator + (const Vector<double_complex>& r1, const Vector<int> &r2);
Vector<double_complex> operator + (const Vector<int>& r1, const Vector<double_complex> &r2);

Vector<double_complex> conj (const Vector<double_complex> &r);

Vector<double> operator * (int x, const Vector<double> &r);
Vector<double> operator * (const Vector<double> &r,int x);
Vector<double> operator * (double x,const Vector<int>& r);
Vector<double> operator * (const Vector<int>& r, double x);

/*Vector<double_complex> operator * (double x,const Vector<double_complex>& r);
Vector<double_complex> operator * (const Vector<double_complex>& r,double x);
Vector<double_complex> operator * (int x,const Vector<double_complex>& r);
Vector<double_complex> operator * (const Vector<double_complex>& r,int x);
*/

Vector<double_complex> operator * (double_complex x,const Vector<double>& r);
Vector<double_complex> operator * (const Vector<double>& r, double_complex x);
Vector<double_complex> operator * (double_complex x,const Vector<int>& r);
Vector<double_complex> operator * (const Vector<int>& r, double_complex x);

Vector<double_complex> operator / (const Vector<double>& r, double_complex x);
Vector<double> operator / (const Vector<double>& r, int x);
Vector<double_complex> operator / (const Vector<double_complex>& r, double x);
Vector<double_complex> operator / (const Vector<double_complex>& r, int x);
Vector<double> operator / (const Vector<int>& r, double x);
Vector<double_complex> operator / (const Vector<int>& r, double_complex x);


double_complex operator * (const Vector<double>& a, const Vector<double_complex> &b);
double_complex operator * (const Vector<double_complex>& a, const Vector<double> &b);
double operator * (const Vector<double>& a, const Vector<int> &b);
double operator * (const Vector<int>& a, const Vector<double> &b);
double_complex operator * (const Vector<double_complex>& a, const Vector<int> &b);
double_complex operator * (const Vector<int>& a, const Vector<double_complex> &b);

Vector<double_complex> operator % (const Vector<double_complex> &a, 
                                   const Vector<double> &b);
Vector<double_complex> operator % (const Vector<double> &a, 
                                   const Vector<double_complex> &b);
Vector<double_complex> operator % (const Vector<double_complex> &a, const Vector<int> &b);
Vector<double_complex> operator % (const Vector<int> &a, const Vector<double_complex> &b);
Vector<double_complex> operator % (const Vector<double> &a, const Vector<int> &b);
Vector<double_complex> operator % (const Vector<int> &a, const Vector<double> &b);


double *conv2double (int Anz, Vector<double_complex> *r);

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
 return (x<0) ? -1 : 1;
}
Vector<double_complex> *conv2vector (int Anz, double *r);
Vector<double_complex> grad2d (double (*f(Vector<double_complex>)),Vector<double_complex> x, double dx);
Vector<double_complex> convd2c (const Vector<double> &r); 
Vector<double> real(Vector <double_complex> r);
Vector<double> imag(Vector <double_complex> r);
double sDreieck (Vector<double>,Vector<double>,Vector<double>);
double sViereck (Vector<double>,Vector<double>,Vector<double>,Vector<double>); 
void dreh (Vector<double> &r, double phi);
Vector<double> dreh (const Vector<double> &r, double dtheta, double dphi);
Vector<double_complex> vdc (const Vector<double> &r);
Vector<double_complex> makeReal (const Vector<double_complex> &r);
Vector<double> cart2sphere(double x, double y, double z);
Vector<double> sphere2cart (double r, double theta, double phi);
void getKSystem (const Vector<double> &n, const Vector<double> &k,
                Vector<double> &e1, Vector<double> &e2, Vector<double> &e3);

Vector<double> arg (Vector<double_complex> &r);
double_complex asin (double_complex z);
double_complex tan  (double_complex z);
double_complex ihoch (int l);
void sphereunitv (Vector<double> &P, Vector<double> &er, Vector<double> &etheta, Vector<double> &ephi); 

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
#define czero Vector<double_complex> (0.0,0.0,0.0)
#define cone Vector<double_complex>  (1.0,1.0,1.0)

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
#endif

#ifndef COMPLEX_I
#define COMPLEX_I
const double_complex I=double_complex (0.0,1.0);
#endif
