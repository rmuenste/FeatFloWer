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
#ifndef KEGELSTUMPF_H
#define KEGELSTUMPF_H

#include "form.h"

/**
	@author Thomas Weigel <weigel@lat.ruhr-uni-bochum.de>
*/
class Kegelstumpf : public Form
{
public:
    Kegelstumpf(const Vector< double >& P, double rmin, double rmax, double h, complex< double > n, 
    double r0=1.0, 
    const Matrix< complex < double > > alpha=CUNITY, 
    const Vector< double >& Ex=ex, 
    const Vector< double >& Ey=ey, 
    const Vector< double >& Ez=ez);
    Kegelstumpf (const Kegelstumpf &E);
    Kegelstumpf(const Form& F);
    Vector<double> norm(const Vector<double> &Ps);
    void binWrite (ofstream &os);
    void binRead (ifstream &is);
    ~Kegelstumpf();
    double schneideMantel (Vector<double> &p, Vector<double> &k);
    bool isInside(const Vector< double >& Ps);
    bool next(const Vector< double >& p, const Vector< double >& k, Vector< double >& pout, const int inside); 
    double schneideDeckel(int i, Vector<double> &p, Vector<double> &k);
    void setParms(double rmin, double rmax, double h)
    {
     if (rmin>rmax) {this->rmax=rmin; this->rmin=rmax;}
     else
     {
     this->rmax=rmax;
     this->rmin=rmin;
     }      
     this->h=h;
     hs=h*this->rmax/(this->rmax-this->rmin);
    }
    void initQuad();
    double getRmax(double rmax) {return rmax;}
    double getRmin(double rmin) {return rmin;}
    void setRmax(double rmax) {this->rmax=rmax;hs=this->rmax/(this->rmax-this->rmin);}
    void setRmin(double rmin) {this->rmin=rmin;hs=this->rmax/(this->rmax-this->rmin);}
    void setr0(double r0);
     void setP (Vector<double> r) {P=r; initQuad();} //
    void setP(double x, double y, double z) {setP(Vector<double> (x,y,z));}     
    void scale (double sf);
    double Volume();
    double rmin; // Kleiner Radius
    double rmax;
    double h;
    double hs; // Höhe des gesamten (nicht abgeschnittenen) Kegels
};

#endif
