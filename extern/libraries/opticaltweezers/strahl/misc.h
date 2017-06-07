
/***************************************************************************
                          misc.h  -  description
                             -------------------
    begin                : Mit Mär 12 2003
    copyright            : (C) 2003 by Thomas Weigel
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
#include "superellipsoid.h"
#include "surface.h"
#include "ellipsoid.h"
#include "funsurf.h"
#include "linse.h"
#include "zylinder.h"
#include "kegel.h"
#include "kegelstumpf.h"
#include "kegelstumpfhohl.h"
#include "compound.h"
#include "superellipsoid_n.h"
#include "erythrocyte.h"
#include "vector.h"
#include "strahl.h"

#ifndef MISC_H
#define MISC_H

#ifndef c_light
#define c_light 299792458.0   // Vakuumlichtgeschwindigkeit in m/s
#endif

#define Z0 376.730313461    // Impedanz des Vakuums Z0=sqrt (mu_0/epsilon_0)=mu_0*c_light in Ohm 
#define mu0 4.0*M_PI*1E-7   // µ0 in N/A^2
#define eps0 8.854187817E-12 // Epsilon_0
#define Planck_h 6.62606896E-34 // Plancksches Wirkungsquantum in Js
#define Planck_hquer 1.054571628E-34 // h/2PI in Js 

using namespace std;

void initInc(Form *E);
void setR0 (Form *E, double r0);
void copyFormList (Form **&d ,Form **s, int anz);
void binWriteIncList (ofstream &os, Form **E, int anz);
void binWriteInc (ofstream &os, Form *E);
void binReadIncList (ifstream &is, Form **&E, int anz);
void binReadInc (ifstream &is, Form *&E, bool isNew);
void copyInc (Form *&d, Form *s);
void deleteInc (Form *E);
Vector<double> force (Vector<double> norm, Ray Se, Ray Sr, Ray St, double df);
double gaussw(double z, double wvl, double w0);
complex<double> gaussphase (Vector<double> P, Vector<double> F, Vector<double> k, double w0, double k0);
double NA2w0(double lambda, double NA, complex<double> n);
float readLE_float32(istream &is);
int readLE_int32(istream &is);
#endif