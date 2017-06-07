include "Detector.h"
#include "matrix.h"




Detector::Detector(void)
{
  n1=0;
  n2=0;
  D=0;
}

Detector::Detector(int n1, int n2)
{
 init(n1,n2);
}

void Detector::init(int n1,int n2)
{
  D=new Vector<complex<double> > *[n1];
 for (int i=0; i<n1; i++)
	 D[i]=new Vector<complex<double> >[n2];
 this->n1=n1;
 this->n2=n2;
}

Detector::Detector (const Detector& Det)
{
 if (Det.n1>0)
 {
	 init(Det.n1,Det.n2);
     for (int i1=0; i1<n1; i1++)
		 for (int i2=0; i2<n2; i2++)
			 this->D[i1][i2]=Det.D[i1][i2];
 }
 n=Det.n;
 e1=Det.e1;
 e2=Det.e2;
 P=Det.P;
}

Detector& Detector::operator = (const Detector& Det)
{
 if (Det.n1>0)
 {
	 init(Det.n1,Det.n2);
     for (int i1=0; i1<n1; i1++)
		 for (int i2=0; i2<n2; i2++)
			 this->D[i1][i2]=Det.D[i1][i2];
 }
n=Det.n;
 e1=Det.e1;
 e2=Det.e2;
 P=Det.P;
 return *this;
}


Detector::~Detector(void)
{
 clear();
}

void Detector::clear()
{
	 if (n1>0) 
  {
	  for (int i1=0; i1<n1; i1++)
		  delete[] D[i1];
      delete[] D;
	  D=0;
	  n1=0;
	  n2=0;
  }
}

int Detector::N1() { return n1; }
int Detector::N2() { return n2; }
void Detector::saveabs (char *fn, int type)
{
  ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=abs(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << endl;
   }
  os.close();
}


void Detector::savereal (char *fn, int type)
{
  ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=real(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << endl;
   }
  os.close();
}


void Detector::saveimag (char *fn, int type)
{
  ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
     if (type>=3) h=abs(D[i1][i2]);
	 else
	  h=imag(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << endl;
   }
  os.close();
}


void Detector::savePhase (char *fn, int type)
{
  ofstream os;
  double h;
  os.open (fn);
   for (int i1=0; i1<n1; i1++)
   {
    for (int i2=0; i2<n2; i2++)
	{
	  h=arg(D[i1][i2][type]);
	 os << h << "   ";
	}
	os << endl;
   }
  os.close();
}

bool Detector::cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l)
{ 
  switch (type)
  {
  case DETECTOR_PLANE : return ((DetectorPlane *)this)->cross(P,k,i1,i2,l); break;
  }
  return false;
}


DetectorPlane::DetectorPlane(void)
{
  n1=0;
  n2=0;
  D=0;
  type=DETECTOR_PLANE;
}

DetectorPlane::DetectorPlane (Vector<double> P, Vector<double> e1, Vector<double> e2, int n1, int n2)
{
    D=new Vector<complex<double> > *[n1];
 for (int i=0; i<n1; i++)
	 D[i]=new Vector<complex<double> >[n2];
 this->n1=n1;
 this->n2=n2;
 this->P=P;
 for (int i=0; i<3; i++)
 {
  if (e1[i]!=0) 
  this->e1[i]=1.0/e1[i]*n1;
  else this->e1[i]=0;
  
  if (e2[i]!=0) 
  this->e2[i]=1.0/e2[i]*n2;
  else this->e2[i]=0;
 }
 this->n=e1%e2;
 this->n=this->n/abs(this->n);
 type=DETECTOR_PLANE;
}

bool DetectorPlane::cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l)
{
 /*Matrix<double> M(e1,e2,-k);
 Vector<double> L;
 bool invertierbar; 
 M=invert(M,invertierbar);
 if (!invertierbar) return false;
 L=M*(P-this->P);
 i1=(int)(L[0]/abs(e1));
 i2=(int)(L[1]/abs(e2));
 Vector<double> h=P-this->P;
 if ( (i1<0) || (i1>=n1) || (i2<0) || (i2>=n2) || (L[2]<0)) return false;
 return true;*/
 double kn=k*n;
 Vector<double> d=this->P-P;
 if (kn==0) return false;
 l=d*n/kn;
 Vector<double> Ph=l*k-d;
 i1=Ph*e1;
 i2=Ph*e2;
  // cout << "i1=" << i1 << "   i2=" << i2 << "   Ph=" << P+l*k << "   l=" << l <<endl;
  if ( (i1<0) || (i1>=n1) || (i2<0) || (i2>=n2)) return false;
 return true;
}

ostream& operator << (ostream &os, Detector& D) 
{
 for (int i1=0; i1<D.n1; i1++)
 { 
  for (int i2=0; i2<D.n2; i2++)
     os << D.D[i1][i2] << "  ";
  os << endl; 
 }
 return os;
}
 

ostream& operator << (ostream &os, mem2D<Vector<complex<double> > > &m)
{
	int n1 = m.getn();
	int n2 = m.getm();
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++) os << m(i, j) << "   ";
		os << endl;
	}
	return os;
}
