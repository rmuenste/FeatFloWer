#ifndef PSTRAHL_H
#define PSTRAHL_H
#include "istrahl.h"
/** Klasse zur Berechnung von Strahlen, bei denen nur die Leistung vorgegeben ist (wichtig z.B. bei den Berechnungen zur optischen Pinzette */
class Ray_pow :
	public IStrahl
{
public:
	Ray_pow(void);
	Ray_pow(double pow, const Vector<double> &p,
         const Vector<complex<double> > &Pol, const Vector<double> &K,
         complex<double>  n0, double r0, double k0,
         const int Anzein, Form **Einschluss);
	Ray_pow(Ray_pow &r) : IStrahl(r)
	{
		this->Pow=r.Pow;
	}
	Ray_pow reflect(Vector<double> n, complex<double> n1, complex<double> n2);
	void refract(Matrix<complex<double> > FT, Vector<double> N, complex<double> n1, complex<double> n2);
	~Ray_pow(void);
	double Pow;
	friend ostream& operator << (ostream &os,Ray_pow S);
};
#endif
