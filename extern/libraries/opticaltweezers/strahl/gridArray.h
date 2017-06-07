#pragma once
#include "form.h"
#include "vector.h"
#include <iostream>
using namespace std;

template <class T> 
class gridArray
{
public:
   gridArray(void);
   gridArray(int ntheta, int nphi, Form *E=0);
   gridArray(gridArray &A);
   gridArray& operator = (gridArray r);
   T operator () (Vector<double> r);
   void clear();
   void init (int ntheta, int nphi);
   ~gridArray(void) {}
inline int getNphi() { return nphi; }
inline int getNtheta() { return ntheta; }
protected: 
	int ntheta,nphi;
	double dphi,dtheta;
	Form *E;
	T **D;	
	bool normal;
};



emplate <class T> gridArray<T>::gridArray(int ntheta, int nphi, Form *E=0)
	{
    D=0;
	normal=false;
	init(ntheta,nphi);
	if (E!=0) normal=true;	
    }


template <class T> gridArray<T>::gridArray(void)
{
	ntheta=0;
	nphi=0;
	D=0;
	E=0;
}



template <class T> void gridArray<T>::clear()
{
  if (D!=0)
  {
	  for (int itheta=ntheta-1; itheta>=0; itheta--)
		  delete[] D[itheta];
	  delete[] D;
	  ntheta=0;
	  nphi=0;
	  D=0;
  }
}

template <class T> gridArray<T>::gridArray(gridArray<T> &A)
{
  D=0;
  init(A.getNtheta(),A.getNphi());	 
  for (int itheta=0; itheta<ntheta; itheta++)
	  for (int iphi=0; iphi<nphi; iphi++)
		  D[itheta][iphi]=A.D[itheta][iphi];
}

template <class T> void gridArray<T>::init (int Ntheta, int Nphi)
{
	if (D!=0) clear();
	this->ntheta=Ntheta;
	this->nphi=Nphi;
	D=new T*  [ntheta];
	for (int itheta=0; itheta<ntheta; itheta++)
		D[itheta]=new T [nphi];
	
}

template <class T> T gridArray<T>::operator () (Vector<double> r)
{
	int itheta, iphi;
	
	Vector<double> rc=cart2sphere(r);
	if (!normal)
	switch (E->Type())
	{
	 case ELLIPSOID : itheta=rc[1]/M_PI*(ntheta-1);		              
		              iphi=rc[2]/(2.0*M_PI)*(nphi-1);
					  break;
	}
	else 
	{
		itheta=rc[1]/M_PI*(ntheta-1);		              
		iphi=rc[2]/(2.0*M_PI)*(nphi-1);
	}
	cout << itheta << "  " << iphi << endl;
	return D[itheta][iphi];
}

template <class T> gridArray<T>& gridArray<T>::operator = (gridArray<T> A)
{
	if (this == &A) return *this;
	init(A.getNtheta(),A.getNphi());
    for (int itheta=0; itheta<ntheta; itheta++)
	  for (int iphi=0; iphi<nphi; iphi++)
		  D[itheta][iphi]=A.D[itheta][iphi]; 
	return *this;
}
