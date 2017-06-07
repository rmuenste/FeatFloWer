/***************************************************************************
                          ebene.h  -  description                              
                             -------------------                                         
    begin                : Sat Feb 19 2000                                           
    copyright            : (C) 2000 by Thomas Weigel                         
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


#ifndef EBENE_H
#define EBENE_H


/**
  *@author Thomas Weigel
  */

#include<iostream>
#include<fstream>
#include "vector.h"
#include "matrix.h"

class Ebene {
public: 
	Ebene();
	Ebene(const Vector<double>&, const Vector<double>&, const Vector<double>&);
	void schneideKugel (Vector<double> O, double r);
	void norm();
	void Normalenform ();
  double Abstand (Vector<double> R);
  void toString(char *S);
  void binWrite (ofstream &os);
  void binRead (ifstream &is);
  void dreh(double dr, double dtheta, double dphi);
	friend ostream & operator << (ostream &os, Ebene E)
  {
   os << " P=" << E.P << endl;
   os << "n=" << E.n << endl;
   os << "e1=" << E.e1 << endl;
   os << "e2=" << E.e2 << endl;
   return os;
  }
	~Ebene();
  Vector<double> P,e1,e2,n;	
};

const Ebene Exz=Ebene(Vector<double>(0.0,0.0,0.0),ex,ez);
const Ebene Exy=Ebene(Vector<double>(0.0,0.0,0.0),ex,ey);
#endif
