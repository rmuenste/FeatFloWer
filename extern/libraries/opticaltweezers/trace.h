
#include "PStrahl.h"
#include "LightSrc.h"
#include "ellipsoid.h"
#include "misc.h"
#include <iostream>
#include <fstream>
using namespace std;

#ifndef TRACE_H
#define TRACE_H

#define NO_SPECIAL_OPTS 0
#define SAVE_FORCE  1
#define SAVE_TORQUE 2

#define MAX_RECURSIONS 10 ///< Maximale Anzahl Rekursionen
#define MAX_REFLEXIONS 5

#define RAY_HITS_OBJECT  1
#define RAY_HITS_NOTHING 0 


/// Ausgabe zu Debugzwecken


typedef struct
{
	Vector<double> P;
	Vector<double> Pold;
	Vector<double> F;
	Vector<double> L;
	Vector<double> k;
} PointInfo;


typedef struct
{
	Vector<double> Pos;
	Vector<double> axis;
	double angle;
} DynCalc;

typedef struct
{
	Vector<double> F, L, P;
} ForceData;


void traceoneRay(Ray_pow &ray, Vector<double> *F, Vector<double> *L, complex<double> ns, int &iRC, int AnzReflex= MAX_REFLEXIONS); ///< einen Strahl verfolgen: ray: Strahl, F,L: Listen mit Kraft (F) und Drehmoment (L), iRC: Counter bisher gemachten Reflexionen, AnzReflex: max. Anzahl Reflexionen
void traceoneRay(Ray_pow &ray, int &anzP, Vector<double> **&P, int &anzF, ForceData *&fd,  complex<double> ns, int reflexCounter, int AnzReflex = 3 , bool drawOutgoing=false);
void trace (int nLS, LightSrc **ls, Vector<double> *F,Vector<double> *L); ///< Alle Strahlen verfolgen, nLS: Anzahl Lichtquellen, lsd
void trace(int nLS, LightSrc **ls, int &anzP, Vector<double> **&P, int &anzF, ForceData *&fd, int AnzReflex = MAX_REFLEXIONS, bool drawOutgoing=false);
void findRootF(Vector<double> &r,int objIndex, int nLS, LightSrc **ls, Vector<double> *F,Vector<double> *L);
void doDynamics (int nLS, LightSrc **ls, int nSteps, double dt, double *rho, DynCalc **&Erg);
#endif
