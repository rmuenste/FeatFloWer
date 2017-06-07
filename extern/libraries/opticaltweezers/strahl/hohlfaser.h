#ifndef HOHLFASER_H
#define HOHLFASER_H
#include "zylinder.h"

class Hohlfaser :
	public Zylinder
{
public:
	Hohlfaser(void);
	~Hohlfaser(void);
    Hohlfaser (const Hohlfaser &E);
    Hohlfaser (
             const Vector<double> &P,
             double ri,
			 double ra,
             double h,
             complex<double>  n,
             double r0=1.0,
             const Matrix<complex<double> > alpha=CUNITY,
             const Vector<double> &Ex=ex,
             const Vector<double> &Ey=ey,
             const Vector<double> &Ez=ez
             );
	bool next(const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside);
	double ri,ra; // Innen- und Aussenradius
	Vector<double> norm(const Vector<double> &Ps);
protected:
	 double schneideDeckel(const Vector<double> &P, const Vector<double> &k, Vector<double> &SP);	 
	 double schneideMantel (Vector<double> &p, Vector<double> &k);
};

#endif