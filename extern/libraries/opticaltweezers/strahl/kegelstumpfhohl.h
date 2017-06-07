/***************************************************************************
 *   Copyright (C) 2005 by Thomas Weigel                                   *
 *   weigel@lat.ruhr-uni-bochum.de                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef KEGELSTUMPFHOHL_H
#define KEGELSTUMPFHOHL_H

#include "form.h"

/**
	@author Thomas Weigel <weigel@lat.ruhr-uni-bochum.de>
*/
class KegelstumpfHohl : public Form
{
public:
    KegelstumpfHohl(const Vector< double >& P, 
    double rimin, double rimax, double ramin, double ramax, 
    double h, complex< double > n, 
    double r0=1.0, 
    const Matrix< complex < double > > alpha=CUNITY, 
    const Vector< double >& Ex=ex, 
    const Vector< double >& Ey=ey, 
    const Vector< double >& Ez=ez);
    double Volume();
    void setParms (double rimin, double rimax, double ramin, double ramax, double h);
    KegelstumpfHohl (const KegelstumpfHohl &E);
    KegelstumpfHohl(const Form& F);
    bool next(const Vector< double >& Ps, const Vector< double >&K, Vector< double >& pout, const int inside);
    double schneideDeckel(int i, Vector<double> &p, Vector<double> &k);
    double schneideMantel (int i, Vector<double> &p, Vector<double> &k);
    Vector<double> norm(const Vector<double> &Ps);
    void binWrite (ofstream &os);
    void binRead (ifstream &is);
    void scale (double sf);
	
	double getRamin() {return ramin; }
	double getRimin() {return rimin; }
	double getRamax() {return ramax; }
	double getRimax() {return rimax; }

	void setRamin(double ramin) {this->ramin=ramin; }
	void setRimin(double rimin) {this->rimin=rimin; }
	void setRamax(double ramax) {this->ramax=ramax; }
	void setRimax(double rimax) {this->rimax=rimax; }

    void setr0(double r0);
    void setP (Vector<double> r) {P=r; initQuad();} //
    void setP(double x, double y, double z) {setP(Vector<double> (x,y,z));}
    bool isInside( const Vector<double> &P) {return true;} // ÄNDERN
    void initQuad(){}
    ~KegelstumpfHohl();
    double rimin,rimax,ramin,ramax;  // rimin, rimax : Radien innen  ramin, ramax: Radien außen
    double h,hsi,hsa;               
};

#endif
