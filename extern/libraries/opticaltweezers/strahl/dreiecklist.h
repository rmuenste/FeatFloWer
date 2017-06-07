#ifndef DREIECKLIST_H
#define DREIECKLIST_H
#include<iostream>
#include "vector.h"
#include "funktion.h"
#include<sstream>

#ifndef SW
#define SW 10.0
#endif

#include "surface.h"


//Vector<double> SurfGit[NTHETA][NPHI];
void SurfGitter(Funktion fu,Vector<double> **SG,int ntheta, int nphi);
void saveGitter(Vector<double> **SG, int ntheta, int nphi);
// void makeSurf(Vector<double> **SG, surface  &Su);
surface makeSurf(Vector<double> **SG,int ntheta, int nphi);
bool nullstelle(Funktion f,double &ns);
#endif

