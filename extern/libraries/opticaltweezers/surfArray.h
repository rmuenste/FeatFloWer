#ifndef SURFARRAY_H
#define SURFARRAY_H

#include"form.h"
#include"vector.h"
#include "ellipsoid.h"
#include "surface.h"
// #include <vector>

/**
  Hier können Daten wie z.B. Kräfte oder Drehimpulse, die  Oberflächen der Objekte wirken, abgespeichert werden.
  Die Einteilung der Surface-Objekte werden direkt übernommen, alle anderen Objekte werden diskretisiert (noch nicht
  implementiert !)
*/

class surfArray
{
public:
	surfArray();
	surfArray(int nObj, Form **Objs);  /// Konstruktor mit Übergabe der Objekte
	surfArray(const surfArray &SA); /// Copy-Konstruktor
	surfArray& operator =(const surfArray& SA); 
	~surfArray();
	Vector<double> operator () (int objIndex, int i) { return data[objIndex][i]; }
	Vector<double> getP(int objIndex, int i) { return P[objIndex][i]; } /// Ort des i-ten Gitters
	void addValue(int objIndex, Form *F, Vector<double> v); /// füge ein beliebiges Objekt hinzu
	void addValue(int objIndex, surface *S, int i, Vector<double> v); /// füge ein Surface-Objekt hinzu
	void clearData(); /// Lösche Gitter 
	void genEllipsoidGitter(int objIndex, FormEllipsoid *E); /// generiere Gitter für Ellipsoid-Objekt (noch nicht fertig)
	int  numItems(int iObj) { return N[iObj]; }
	int numObjs() { return nObjs; }
	
	
private:
	Vector<double> **data;
	Vector<double> **P;
		int nObjs;
		int *N;
		int *type;
};

#endif
