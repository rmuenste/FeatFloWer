#ifndef MATRIX_H
#define MATRIX_H
#include "vector.h"
#include <iostream>
/** Template-Klasse, die eine 3x3 Matrix repräsentiert */
template<class T> class Matrix 
{
 public :
 Matrix ()
 {  
 for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
   M[i][j]=0.0; 
 }

 void binWrite (ofstream &os) /// speichert die Matrix binär ab
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
  os.write ((char *) &M[i][j], sizeof (M[i][j])); 
 }

 void binRead (ifstream &is)
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
  is.read ((char *) &M[i][j], sizeof (M[i][j])); 
 }
 
 Matrix (Vector<T> a, Vector<T> b, Vector<T> c)  /// Konstruktor von Matrix, die Vektoren a,b und c stellen die Spalten dar
 {
  for (int i=0; i<3; i++)
     M[i][0]=a[i];
  for (int i=0; i<3; i++)
     M[i][1]=b[i]; 
  for (int i=0; i<3; i++)
     M[i][2]=c[i]; 
 }
 
 Matrix& operator = (const Matrix& A) 
 {
 if (this==&A) return *this;
 for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) 
    M[i][j]=A.M[i][j];
 return *this;
 }

 bool operator == (const Matrix& A)
 {
  bool Erg=true;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
      Erg=Erg && (M[i][j]==A.M[i][j]);
   return Erg;
 } 
 
 bool operator != (const Matrix& M)
 {
  bool Erg=false;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
      Erg=Erg || (M[i][j]!=M.M[i][j]);
   return Erg;
 } 
  
 T& operator () (int i, int j) { return M[i][j]; } 
 
 Matrix operator + (const Matrix& A)
 {
 Matrix<T> Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=M[i][j]+A.M[i][j];
  return Erg;
 }

 Matrix& operator +=(const Matrix& A)
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    M[i][j]+=A.M[i][j];
  return *this;
 }

Matrix operator - ()
 { 
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=-M[i][j];
  return Erg;  
 }

 Matrix& operator -=(const Matrix& A)
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    M[i][j]-=A.M[i][j];
  return *this;
 }
 
 Matrix& operator *= (const Matrix& A)
 {
  T m[3][3];
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
   {
     m[i][j]=0;
     for (int l=0; l<3; l++)
     m[i][j]+=M[i][l]*A.M[l][j];
    }
   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
      M[i][j]=m[i][j];
   return *this;
 }

 Matrix& operator /=(const T& x)
 {
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    M[i][j]/=x; 
   return *this; 
 }

 Matrix& operator *=(const T& x)
 {
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    M[i][j]*=x; 
   return *this; 
 }
 
 
 friend Matrix operator - (const Matrix& A, const Matrix& B)
 {
 Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=A.M[i][j]-B.M[i][j];
  return Erg;
 }
 
 friend Matrix operator * (const Matrix& A, const Matrix& B)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
     for (int l=0; l<3; l++)
      Erg.M[i][j]+=A.M[i][l]*B.M[l][j];
  return Erg; 
 }
 
 friend Matrix operator * (const Matrix& A, const T& x)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=A.M[i][j]*x;
  return Erg; 
 }
 
 friend Matrix operator * (const T& x, const Matrix& A)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=x*A.M[i][j];
  return Erg; 
 }
  
 friend Vector<T> operator * (const Matrix& A, const Vector<T>& r)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    Erg[i]+=A.M[i][j]*r[j];
  return Erg;  
 }
 
 friend Vector<T> operator * (const Vector<T>& r, const Matrix& A)
 {
  Vector<T> Erg;
  for (int i=0; i<3; i++) 
   for (int j=0; j<3; j++)
    Erg[i]+=A.M[i][j]*r[j];
  return Erg;  
 }
 
/* friend Matrix operator / (const Matrix& A, const T& x)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
     for (int l=0; l<3; l++)
      Erg.M[i][j]+=A.M[i][l] / x;
  return Erg; 
 }*/
 
friend Matrix operator / (const Matrix& A, const T& x)
 {
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
      Erg.M[i][j]=A.M[i][j] / x;
  return Erg; 
 }

 friend ostream& operator << (ostream& os, const Matrix& A)
 {
  for (int i=0; i<3; i++)
   os << A.M[i][0] << "  " << A.M[i][1] << "  " << A.M[i][2] << endl;
  return os; 
 }  

 friend istream& operator >> (istream& is, Matrix& A)
 {
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    is >> A.M[i][j];
  return is; 
 }
 
 friend Matrix transpose (const Matrix& A)
 { 
  Matrix Erg;
  for (int i=0; i<3; i++)
   for (int j=0; j<3; j++)
    Erg.M[i][j]=A.M[j][i];   
  return Erg;
 }  
  
 friend T trace(const Matrix& A)
 {
  T Erg=0;
  for (int i=0; i<3; i++)
   Erg+=A.M[i][i];
  return Erg;
 }
 
 friend T det(const Matrix& A)
 {
  T Erg; 
  Erg=A.M[0][0]*(A.M[1][1]*A.M[2][2]-A.M[2][1]*A.M[1][2])
     -A.M[0][1]*(A.M[1][0]*A.M[2][2]-A.M[2][0]*A.M[1][2])
     +A.M[0][2]*(A.M[1][0]*A.M[2][1]-A.M[2][0]*A.M[1][1]);
  return Erg;
 }

 friend Matrix streich (const Matrix& M,int i, int j)
 {
  Matrix Erg;
  for (int k=0; k<3; k++)
   for (int l=0; l<3; l++)
    {
     if ((k!=i)&&(l!=j)) Erg.M[k][l]=M.M[k][l];
     if ((l==j)&&(k==i)) Erg.M[k][l]=1; 
     else
     {
       if (l==j) Erg.M[k][l]=0;
       if (k==i) Erg.M[k][l]=0;
     }
    }
  return Erg; 
 }
 
 friend Matrix invert (const Matrix<T>& A,bool &invertierbar)
 {
  Matrix Erg;
  T C;
  C=det(A);
  if (C==0.0) invertierbar=false;
  else
  {
   Erg.M[0][0]=A.M[1][1]*A.M[2][2] - A.M[2][1]*A.M[1][2];
   Erg.M[1][0]=A.M[2][0]*A.M[1][2] - A.M[1][0]*A.M[2][2];
   Erg.M[2][0]=A.M[1][0]*A.M[2][1] - A.M[2][0]*A.M[1][1];

   Erg.M[0][1]=A.M[2][1]*A.M[0][2] - A.M[0][1]*A.M[2][2];
   Erg.M[1][1]=A.M[0][0]*A.M[2][2] - A.M[2][0]*A.M[0][2];
   Erg.M[2][1]=A.M[2][0]*A.M[0][1] - A.M[0][0]*A.M[2][1];

   Erg.M[0][2]=A.M[0][1]*A.M[1][2] - A.M[1][1]*A.M[0][2];
   Erg.M[1][2]=A.M[1][0]*A.M[0][2] - A.M[0][0]*A.M[1][2];
   Erg.M[2][2]=A.M[0][0]*A.M[1][1] - A.M[0][1]*A.M[1][0];
   invertierbar=true;
  }
 return Erg/C; 
 }

 friend Matrix invert (const Matrix<T>& A)
 {
  bool dummy;
  return invert (A,dummy); 
 } 

 void clear()
 {
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      M[i][j]=0;
 } 

 T M[3][3]; 
};

// Matrix-Matrix Operatoren mit gemischten Typen (complex<double> /double)
Matrix<complex<double> > operator * (const Matrix<double>&, const Matrix<complex<double> >&);
Matrix<complex<double> > operator * (const Matrix<complex<double> >&, const Matrix<double>&);
Matrix<complex<double> > operator + (const Matrix<double>&, const Matrix<complex<double> >&);
Matrix<complex<double> > operator + (const Matrix<complex<double> >&, const Matrix<double>&);
Matrix<complex<double> > operator - (const Matrix<double>&, const Matrix<complex<double> >&);
Matrix<complex<double> > operator - (const Matrix<complex<double> >&, const Matrix<double>&);

// Matrix-Vektor Operatoren mit gemischten Typen (complex<double> /double)
Vector<complex<double> > operator * (const Matrix<double>&, const Vector<complex<double> >&);
Vector<complex<double> > operator * (const Matrix<complex<double> >&, const Vector<double>&);

// Matrix-Skalar Operatoren mit gemischten Typen (complex<double> /double)
Matrix<complex<double> > operator * (const Matrix<double>&, const complex<double> &);
Matrix<complex<double> > operator * (const Matrix<complex<double> >&, const double&);
Matrix<complex<double> > operator * (const double&, const Matrix<complex<double> >&);
Matrix<complex<double> > operator * (const complex<double> &, const Matrix<double>&);
Matrix<complex<double> > operator / (const Matrix<double>&, const complex<double> &);
Matrix<complex<double> > operator / (const Matrix<complex<double> >&, const double&);


// Drehmatrizen um die 3 Raumachsen
Matrix<double> Dx(double phi);
Matrix<double> Dy(double phi);
Matrix<double> Dz(double phi);

Matrix<double> unity();  // Einheitsmatrix
Matrix<complex<double> > cunity (); 
Matrix<double> null(); // Nullmatrix 
Matrix<double> drehmatrix(const Vector<double> a, double gamma); /// Drehmatrix für Drehung um die Achse a um den Winkel gamma
Matrix<double> drehmatrix(double alpha, double beta, double gamma);
Matrix <double> drehmatrixD (Vector<double> n, Vector<double> k, double gamma);
Matrix<double> rotMat(Vector<double> P, double dtheta, double dphi);
void trafo (const Vector<double> &, 
            const Vector<double> &, 
            const Vector<double> &, 
            Matrix<double> &,
            Matrix<double> &);

const Matrix<double> UNITY=Matrix<double>(Vector<double>(1,0,0),Vector<double>(0,1,0),Vector<double>(0,0,1));
const Matrix<complex<double> > CUNITY=Matrix<complex<double> >(Vector<complex<double> >(1,0,0),Vector<complex<double> >(0,1,0),Vector<complex<double> >(0,0,1));

#endif // MATRIX_H
