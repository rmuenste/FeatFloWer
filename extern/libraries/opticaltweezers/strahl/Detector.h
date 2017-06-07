#ifndef DETECTOR_H
#define DETECTOR_H

#include "vector.h"
#define DETECTOR_ANGLE 1
#define DETECTOR_PLANE 2
#include <iostream>
using namespace std;
/*!
  Klasse Detector beschreibt die Grundzüge eines abstrakten "Detektors"
*/

#define MEM2D_DIR_X 0
#define MEM2D_DIR_Y 1
#define MEM2D_DIR_Z 2


template <class T> class mem2D
{
public:
	~mem2D() { clear(); }
	mem2D() { n = 0; m = 0; data = 0; }
	mem2D(int n, int m, double w1, double w2, int dir1=MEM2D_DIR_X, int dir2=MEM2D_DIR_Z) 
	/**
	Initialization 
	n,m : Number of elements in e1- and e2-direction
	w1,w2 : width in e1- and e2-direction
	*/
	{ 
		Vector<double> h1, h2, e1, e2;
		this->n = n;
		this->m = m;
		l1 = w1 / ((double)n - 1);
		l2 = w2 / ((double)m - 1);
		d = Vector<double>(1.0,1.0,1.0);
		switch (dir1)
		{
		case MEM2D_DIR_X: d[0] = l1; e1 = ex;  break;
		case MEM2D_DIR_Y: d[1] = l1; e1 = ey; break;
		case MEM2D_DIR_Z: d[2] = l1; e1 = ez; break;
		}

		switch (dir2)
		{
		case MEM2D_DIR_X: d[0] = l2; e2 = ex; break;
		case MEM2D_DIR_Y: d[1] = l2; e2 = ey; break;
		case MEM2D_DIR_Z: d[2] = l2; e2 = ez; break;
		}

		data = new T*[n];
		for (int i = 0; i < n; i++)
			data[i] = new T[m];
		h1 = e1 / abs(e1);
		h2 = e2 / abs(e2);
		
		d1 = h1/l1;
		d2 = h2/l2;

		P = -w1 / 2.0*e1 - w2 / 2.0*e2; /* Origin (left lower corner) */		
	}
	
	inline mem2D(mem2D& d)
	{
//		clear();
		n = d.getn();
		m = d.getm();
		data = new T*[n];
		for (int i = 0; i < n; i++)
		{
			data[i] = new T[m];
			for (int j = 0; j < m; j++)
				data[i][j] = d(i, j);
		}
	}

	T& operator () (Vector<double> P)
	{
	//	P-P0=i*d1+j*d2;
    // returns data corresponding to point P
		Vector<double> h = (P - this->P);
		int i = h * d1;
		int j = h * d2;
		return data[i][j];
	}

	T& operator () (int i, int j) { return data[i][j]; }

	void clear() 
	{
		if (data!=0)
		for (int j = 0; j < m; j++) delete[] data[j];
		delete[] data;
		n = 0;
		m = 0;
		data = 0;
	}


	Vector<double> next(Vector<double> p0, Vector<double> k0, double eps=1E-10)
	{
		double lambdax, lambday, lambdaz, lambda;
		double signx, signy, signz;
		double sx, sy, sz;

		signx = copysign(1.0, k0[0]); signy = copysign(1.0, k0[1]); signz = copysign(1.0, k0[2]);

		sx = (floor((p0[0] + signx * 2 * eps) / d[0]) + signx + (signx<0)*(fmod(p0[0] + signx * 2 * eps, d[0]) != 0))*d[0];
		sy = (floor((p0[1] + signy * 2 * eps) / d[1]) + signy + (signy<0)*(fmod(p0[1] + signy * 2 * eps, d[1]) != 0))*d[1];
		sz = (floor((p0[2] + signz * 2 * eps) / d[2]) + signz + (signz<0)*(fmod(p0[2] + signz * 2 * eps, d[2]) != 0))*d[2];

		lambdax = (sx - p0[0]) / k0[0];
		lambday = (sy - p0[1]) / k0[1];
		lambdaz = (sz - p0[2]) / k0[2];

		if (lambdax <= lambday)
			lambda = lambdax;
		else
			lambda = lambday;

		if (lambda >= lambdaz)
			lambda = lambdaz;


		return  p0 + lambda*k0;

	}


	/*
	Vector<double> next(Vector<double> P, Vector<double> k, double eps=1E-10)
	{
		Vector<double> Pn;
		double l[4];
		double min;
	// crossing between line r=P+l*k and rows/cols 
    // at first, find the corresponding rows and colsl[0]
		Vector<double> h = (P - this->P);
		int i = h *d1;
		int j = h *d2;
   
   // row: P=P0+j*d2+l1*d1 => P0+j*d2+l1*d1=P+l*k => (*d2) => P0*d2+j=P*d2+l*(k*d2)
   //      ((P0-P)*d2+j)/(k*d2);
		l[0] = ((this->P - P)*d2 + j) / (k*d2);
		l[1] = ((this->P - P)*d2 + (j+1)) / (k*d2);
		l[2] = ((this->P - P)*d1 + i) / (k*d1);
		l[3] = ((this->P - P)*d1 + (i + 1)) / (k*d1);

		if (l[0] > eps) min = l[0];
		if ((l[1] < min) && (l[1]>eps) ) min = l[1];
		if ((l[2] < min) && (l[2]>eps) ) min = l[2];
		if ((l[3] < min) && (l[3]>eps) ) min = l[3];
		cout << "P=" << P << "   min=" << min << endl;
		return P + min*k;
	}
	*/
	

	void storeE(Vector <double> from, Vector<double> to, Vector<complex<double> >E, Vector<double> k, complex <double> n, double wvl)
	{
		int i, j;
		double l,d,dl;
		Vector<double> h,P,Pn;
		Vector<complex<double> > Ef;
		double k0 = 2.0*M_PI / wvl;
		d = abs(to - from);
		l = 0;
		Ef = E;
		do
		{
			P = from + l*k;
			Pn = next(P, k);
		//	cout << P << "  " << Pn << endl;
 			dl = abs(Pn - P);
			Ef = Ef*exp(I*k0*dl*n);
			l += dl;			
			h = (P - this->P);
			i = h * d1;
			j = h * d2;
		//	cout << "i=" << i << "   j=" << j << "   Ef=" << Ef << endl;
			data[i][j]=data[i][j]+Ef;
		} while (l < d);
	}

	
	
	mem2D& operator = (const mem2D& d)
	{
		if (this == &d) return *this;
		clear();
		n = d.getn();
		m = d.getm();
		data = new T[n];
		for (int i = 0; i < n; i++)
		{
			data[i] = new T[m];
			for (int j = 0; j < m; j++)
				data[i][j] = d(i, j);
		}
		return *this;
	}

	
	int getn() { return n; }
	int getm() { return m; }
	Vector<double> getP() { return P; }
	Vector<double> getd1() { return d1; }
	Vector<double> getd2() { return d2; }
	double getw1() { return w1; }
	double getw2() { return w2; }
	 
 protected: 
	 T **data; // Array for data
	 int n, m; // No. of rows and cols
	 Vector<double> P; // Vector to the lower left corner
	 Vector<double> d1, d2; /// Vectors for one step in cols and rows  (orthogonality of d1 and d2 is assumed !)
	 double l1, l2; /// length of one cols/rows
	 Vector<double> d; /// Vector with the grid width in the corresponding direction (x,y or z)
};

class Detector 
{
public:
	Detector(void);
	Detector(int n1, int n2); /// Detektor mit einem Gitter n1 x n2
	Detector(const Detector &Det);
	Detector& operator= (const Detector& D);
	Vector<complex<double> >& operator () (int i1,int i2) { return D[i1][i2]; }
	~Detector(void);
	void init (int n1, int n2); ///< initialisiere Datenmatrix D
	void clear(); ///< lösche Daten
	bool cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l);
	int N1();
	int N2();
	int Type() {return type; } ///< Art des Detektors
	void saveabs (char *fn, int type=3);  ///< zum Abspeichern des Betrags von E (type=3) oder einer seiner kartesischen Komponenten (type: x->0, y->1, z->2)
	void savePhase (char *fn, int type); ///< zum Abspeichern der Phaseninformation 
    void savereal (char *fn, int type); ///< zum Abspeichern des Realteils der type-ten Komponente
    void saveimag (char *fn, int type); 


friend ostream& operator << (ostream &os, Detector& D);  
Vector<double> e1, e2, P,n;    
	protected: 
Vector<complex<double> > **D; ///< Hier werden die Daten gespeichert (elektrisches Feld)
int n1,n2;
int type;
friend class DetectorAng; 
friend class DetectorPlane;

}; 

class DetectorAng : public Detector
{
};

class DetectorPlane : public Detector
{
public:
  DetectorPlane (void);
  DetectorPlane (Vector<double> P, Vector<double> e1, Vector<double> e2, int n1, int n2);
  bool cross(Vector<double> P, Vector<double> k, int &i1, int &i2, double &l);	
};

#define SAVE_X 0
#define SAVE_Y 1
#define SAVE_Z 2
#define SAVE_PHASE_X 3
#define SAVE_PHASE_Y 4
#define SAVE_PHASE_Z 5
#define SAVE_ABS 6



void save(mem2D<Vector<complex<double> > > m, char *fname, int type = SAVE_PHASE_Y)
{
	int n1 = m.getn();
	int n2 = m.getm();
	ofstream os;
	os.open(fname);
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
			switch (type)
		{
			case SAVE_ABS: os << abs(m(i, j)) << "  "; break;
			case SAVE_PHASE_X: os << arg(m(i, j)[0]) << "  "; break;
			case SAVE_PHASE_Y: os << arg(m(i, j)[1]) << "  "; break;
			case SAVE_PHASE_Z: os << arg(m(i, j)[2]) << "  "; break;
			case SAVE_X: os << abs(m(i, j)[0]) << "  "; break;
			case SAVE_Y: os << abs(m(i, j)[1]) << "  "; break;
			case SAVE_Z: os << abs(m(i, j)[2]) << "  "; break;
		}

		os << endl;
	}
	os.close();
}

ostream& operator << (ostream &is, mem2D<Vector<complex<double> > > &m);
#endif