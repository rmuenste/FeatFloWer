#ifndef ASDD 
#define ASDD 
#include "matrix.h"
#include "vector.h"
 /** 
   Berechnet den Schnbittpunkt einer Geraden mit einem Rechteck im dreidimensionalen Raum 
   Rechteck: Pr: Ecke, e1,e2: Vektoren in Richtung der Kanten und mit deren Länge  r=Pr+L[0]*e1+L[1]*e2
   Gerade :  Ps: Aufpunkt, k: Richtungsvektor       r=Ps+L[2]*k 
*/
  Vector<double> intersectLineRect(Vector<double> Pr,Vector<double> e1, Vector<double> e2, Vector<double> Ps, Vector<double> k, Vector<double> &L);

#endif

