#include "LightSrc.h"
#include "linse.h"
#include "ellipsoid.h"
#include "surface.h"
#include "trace.h"
#include <fstream>
#include <iostream>
using namespace std;
/**************************************************************************************

  Beispielprogramm zur Kraftberechnung:
  Kugel mit Radius 30 wird durch einen fokussierten Gausstrahl hindurchbewegt
  dabei werden die Kraefte und Drehmomente berechnet 

 **************************************************************************************/ 



int main(int argc, char **argv)
{
  Form **O = new Form*[1];          // Da stehen die Objekte drin (hier nur eins)
  LightSrc **L = new LightSrc*[1];  // Da stehen die Lichtquellen drin
  int N = 200;                      // Anzahl Strahlen  
  double r0 = 20000;                // Groesse der "Weltkugel" (=Berechnungsbereich)
  Vector<double> StartPos(0,0,-100);  // Position der Lichtquelle
  double wvl=1.0;                   // Wellenlaenge
  double w0=0.01;                   // (fiktiver) Strahldurchmesser im Fokus  
  Vector<double> focuspos=zero;     // Fokusposition (zero = Nullvektor)     
  double R=30;                      // Radius des Objekts  
  complex<double> n=1.5;            // Brechungsindex des Objekts   

  O[0] = new FormEllipsoid(zero, Vector<double>(R, R, R), n, r0); // Objekt (Ellipsoid) erstellen

  L[0] = new LightSrcGauss(StartPos, N, wvl, w0, focuspos, 100.0); // Lichtquelle erstellen (Gausstrahl)
  L[0]->setN0(1.0); // Brechungsindex Umgebungsmedium setzen
  L[0]->ObjectList(1, O); // Objektliste an Lichtquelle uebergeben
  ((LightSrcGauss*)L[0])->setNA(0.99); // numerische Apertur des Strahls ändern (=sin (Oeffnungswinkel))
  cout << "w0=" <<    ((LightSrcGauss*)L[0])->w0 << endl; // Durchmesser im Fokus hat sich geändert

  Vector<double> **P; 
  Vector<double> *F=new Vector<double> [1]; // Liste mit Kräften (1 wegen nur einem Objekt)
  Vector<double> *D=new Vector<double> [1]; // Liste mit Drehmomenten ( " )

  ofstream os;
  os.open ("test.dat"); 
  for (double x=-1.5*R; x<1.5*R; x+=3.0*R/101.0) 
  {
    L[0]->reset(); // Lichtquelle zuruecksetzen
    L[0]->Ein[0]->P=x*ex;
    F[0]=zero; 
    D[0]=zero;   
    trace (1,L,F,D); // eigentliche Strahlverfolgung
    cout  << L[0]->Ein[0]->P << "   " << F[0] << "   " << D[0] << endl;
    os << x << "   " << F[0] << "   " << D[0] << endl;
  }
  os.close();

  delete[] F;
  delete[] L; 
  delete[] O;
  delete[] D;
}
