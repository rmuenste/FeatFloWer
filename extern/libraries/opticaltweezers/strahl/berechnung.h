#include <complex>
#include <math.h>
#include "vector.h"
#ifdef WITH_SUPERGITTER
#include "supergitter.h"
#endif
#include "gitter.h"
#include "ebene.h"

/** Berechnet den Punkt ix,iy eines Einstrahlgitters, welches sich im Abstand r
    unter einem Winkel theta,phi vor dem Partikel befindet
*/
Vector<double> startpunkt(int ix, int iy,  gitter& git, 
                          double r, double theta, double phi);

/** Berechnet den Punkt ix,iy eines Einstrahlgitters, welches sich
    unter einem Winkel theta,phi direkt vor dem Partikel befindet
*/
Vector<double> startpunkt(int ix, int iy, gitter& git,
                          double theta, double phi);

/// Dreht den Vektor vein um theta und um phi
Vector<double> drehvektor(Vector<double> vein, double theta, double phi);

/// Dreht den Vektor vein um phi
Vector<double> drehphivektor(Vector<double> vein, double phi);

/** Gibt den Ursprung eines Einstrahlgitters wieder, dessen Mittelpunkt um  
    Winkel theta,phi gedreht wurde und das ursprünglich um dr über der Ebene 
    z = git.zmax lag. Wird in der Routine startpunkt benoetigt
*/
Vector<double> ursprungrot(double dr, double theta, double phi, 
                           Vector<double> &exrot, Vector<double> &eyrot, 
                           gitter &git);

/** Berechnet die eigentliche Strahlverfolgung im Gitter git bei Einstrahlung  
 an dem Startpunkt p0  in Richtung k0
 */
#ifdef WITH_SUPERGITTER
SuperGitter verfolgung(Vector<double> p0, Vector<double> k0, SuperGitter &git);


Vector<double> pnext(Vector<double> p0, Vector<double> k0, SuperGitter &git);

Vector<double> pnext(Vector<double> p0, Vector<double> k0, SuperGitter &git, double eps);
#endif
Vector<double> pnext(Vector<double> p0, Vector<double> k0,
Vector<double> d, double eps);
Vector<double> pnext(Ebene E, Vector<double> p0s, Vector<double> k0s,
Vector<double> d, double eps);

