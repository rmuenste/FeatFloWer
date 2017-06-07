#include "vector.h"

class Spherical : public Vector<double>
{
 public :

 Spherical (); 
 Spherical (double,double,double);
 Spherical (const Vector<double>&); 
 Spherical (double,double); 
 Spherical& operator = (const Vector<double> &);
 Vector<double> er();
 Vector<double> ephi();
 Vector<double> etheta();
};
