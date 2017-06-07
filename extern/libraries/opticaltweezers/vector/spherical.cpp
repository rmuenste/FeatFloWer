#include "spherical.h"

Spherical::Spherical ()
{
 data[0]=1.0;
 data[1]=0.0;
 data[2]=0.0;
}

Spherical::Spherical (const Vector<double>& r)
{
 Vector<double> s; 
 s=cart2sphere (r);
 data[1]=s[1];
 data[2]=s[2];
 data[0]=1.0;
}

Spherical::Spherical (double x, double y, double z)
{
 Vector<double> s; 
 s=cart2sphere (Vector<double>(x,y,z));
 data[1]=s[1];
 data[2]=s[2];
 data[0]=1.0;
}
 
Spherical::Spherical(double theta, double phi)
{
 data[0]=1.0;
 data[1]=theta;
 data[2]=phi;
}

Spherical& Spherical::operator = (const Vector<double>& S)
{
 if (&S==this) return *this;
 for (int i=0; i<3; i++)
  data[i]=S[i]; 
 return *this; 
}

Vector<double> Spherical::er()  
{
 Vector<double> Erg;
 Erg[0]=sin(data[1])*cos(data[2]);
 Erg[1]=sin(data[1])*sin(data[2]);
 Erg[2]=cos(data[1]);
 return Erg;
}

Vector<double> Spherical::etheta()
{
 Vector<double> Erg;
 Erg[0]=cos(data[1])*cos(data[2]);
 Erg[1]=cos(data[1])*sin(data[2]);
 Erg[2]=-sin(data[1]);
 return Erg;
}

Vector<double> Spherical::ephi()
{
 Vector<double> Erg;
 Erg[0]=-sin(data[2]);
 Erg[1]=cos(data[2]);
 Erg[2]=0.0;
 return Erg;
}
