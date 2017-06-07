#include "intersect.h"

Vector<double> intersectLineRect (Vector<double> Pr,Vector<double> e1, Vector<double> e2, Vector<double> Ps, Vector<double> k, Vector<double> &L)
{ 
 Matrix<double> M(e1,e2,-k);
 bool inv;
 Matrix<double>  Mi=invert(M,inv);
 if (inv)
 {
 Vector<double> P=Ps-Pr;
 L=Mi*P;
//  cout << "intersect:::L=" << L << "   Ps=" << Ps << "   k=" << k << "   Pr=" << Pr << "   e1=" << e1 << "   e2="  << e2 << endl;

 if ((L[0]<0.0) || (L[0]>1.0) || (L[1]<0.0) || (L[1]>1.0)) return nanV("");
 return Ps+L[2]*k;
 }
}

Vector<double> intersectLineLine( Vector<double> P1, Vector<double> k1, Vector<double> P2, Vector<double> k2)
{

  Vector<double> t;
  return t;
 
}

