#include "vector.h"
#include "zylinder.h"
#include "istrahl.h"
#include <iostream>
#include <fstream>
using namespace std;

int main (int argc, char **argv)
{
 ofstream os;
 os.open ("test.dat");
 ZylinderHexagonal Z(10,20); 

 cout << "Z=" << Z << endl;
 Vector<double> P,Ps,k;
 IStrahl is(Vector<double>(0.0,0.0,0.0),Vector<complex<double> >(0.0,0.0,1.0),Vector<double> (0.0,1.0,0.0),1.0,10000.0,1000.0);
  int i;
  double l;
k=Vector<double>(0,-1,0);
is.P=Vector<double> (0.0,300.0,0.0);
l=Z.schneideMantel(P,k);
cout << "l=" << l  << endl;
for (double x=-10; x<=10; x+=0.01)
 for (double z=0; z<=20; z+=0.5)
  {
   P=Vector<double> (x,200.0,z);
   is.P=P;
   is.k=k;
//   Ps=is.intersectRect(Vector<double> (10.0,0.0,20.0),20.0*ex,40.0*ez);
   l=Z.schneideMantel(P,k);

/*   os << P << "   "  << Ps << endl;
   cout << P << "    " << Ps << endl;*/
 os << P << "   " << P+l*k << endl;
   cout << P << "   " << P+l*k << endl; 
 }



/*
 while (true)
  { 
   for (i=0; i<3; i++)
   {
    cout << "P[" << i << "]=";
    cin >> P[i];
   }
   
  for (i=0; i<3; i++)
   {
    cout << "k[" << i << "]=";
    cin >> k[i];
   }
   if (k[0]<0) break;
   l=Z.schneideMantel(P,k);
   os << P << "   " << P+l*k << endl;
   cout << P << "   " << P+l*k << endl;

}*/
 os.close();
 return 0;
}
