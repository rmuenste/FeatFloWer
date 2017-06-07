#ifndef ERYTHROCYTE_H
#define ERYTHROCYTE_H

#include <lgsurface.h>

/**
	@author Thomas Weigel <weigel@lat.rub.de>
*/

class Erythrocyte : public LGSurface
{
public:
    Erythrocyte();
    Erythrocyte(const Erythrocyte &E);
     double f(const Vector<double> &P) ;     
    Vector<double> grad(const Vector<double> &P);
    double calcG(double t1, double t2);
    ~Erythrocyte();
    void binRead (ifstream &is);
    void binWrite (ofstream &os);
double d,a,b,c,h;
protected :
void init();
};

#endif
