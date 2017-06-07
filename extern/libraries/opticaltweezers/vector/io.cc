#include <iostream>
#include <fstream>
#include <math.h>
#include "resutil.h"


void read_input (char *FName, GlobalParms &Parms)
{
 // liest aus der Datei FName die globalen Parameter Parms aus
 /*
    vergleiche Definition von GlobalParms in >resutil.h<
 */

 ifstream is;
 int dummy;
  
 is.open(FName);
 //Parms.l0=318.4E-9; 
 Parms.l0=2.0*M_PI;  
 is >> Parms.nx >> Parms.ny;    
 is >> Parms.pol >> dummy;      
                                // dummy=1 -> Ausgabedatei enthaelt Phasenbild 
 Parms.phase=(dummy==1); 
 
 is >> Parms.AnzReflex >> Parms.AnzRays; 
 is >> Parms.ResRad >> Parms.ResAzi;      
 is >> Parms.n0;
 is >> Parms.AngleTol >> Parms.evan;
 //Parms.r0=Parms.l0*GeoX(Parms.pol,Parms.ResRad,Parms.n0,1.0,Parms.ResAzi)/(2.0*M_PI);
 Parms.alpha=0.0;
 //Parms.bmax=Parms.r0*Parms.n0; 
 Parms.k0=2.0*M_PI/Parms.l0; 
 Parms.dy=2.0*Parms.r0/Parms.ny;
 Parms.dx=2.0*Parms.r0/Parms.nx;
 Parms.dxy=sqrt(Parms.dx*Parms.dx+Parms.dy*Parms.dy);
 //Parms.db=Parms.n0*Parms.r0/Parms.AnzRays;
 //Parms.phi=2.0*asin(Parms.dx*Parms.dy/(2.0*Parms.r0)); 
 cout << "# Parameter :" << endl;
 cout << "# Resonanz  :";
 if (Parms.pol==SENKRECHT) cout << "TE  ";
 else cout << "# TM  ";  
 cout << Parms.ResRad << "  " << Parms.ResAzi << endl;
 cout << "# Radius r0=" << Parms.r0 << "  Wellenl. l0=" << Parms.l0 << endl;
 cout << "# Miep.  x0=" << 2.0*M_PI*Parms.r0/Parms.l0 << "  k0=" << Parms.k0 << endl;
 cout << "# Anz. Reflexionen: " << Parms.AnzReflex << endl;
 cout << "# " << Parms.nx << " x " << Parms.ny << " Punkte" << endl; 
 cout << "# Anz. Strahlen   :" << Parms.AnzRays << endl;
 cout << "# Brechungsindex  :" << Parms.n0 << endl;
 cout << "# Polarisation    :";  
 if (Parms.pol==SENKRECHT) 
   cout << "  SENKRECHT" << endl;
 else
   cout << "  PARALLEL" << endl;
 cout << "# Winkel-Toleranz :" << Parms.AngleTol << endl;
}

void read_inc (char *FName, int &Anzahl, EinschlussInfo *Ein, double r0)
{
 /*
    liest aus der Datei >FName< die Informationen ueber die 
    Einschluesse aus 
    (vgl. Def. von EinschlussInfo in >resutil.h<
 */  

 ifstream is;
 is.open(FName);
 is >> Anzahl;  
 Ein=new EinschlussInfo[Anzahl+2];
 for (int i=0; i<Anzahl; i++)
 {
  is >> Ein[i].P[0] >> Ein[i].P[1] >> Ein[i].P[2] >> Ein[i].a >> Ein[i].n; 
  cout << "# Eins. " << i << ": ( " 
       << Ein[i].P[0] << "," 
       << Ein[i].P[1] << ","
       << Ein[i].P[2] << " )  " 
       << Ein[i].a << "  " << Ein[i].n << endl; 
  Ein[i].P*=r0;
  Ein[i].a*=r0; 
 }
}


void info (void)
{
 cout << "Benutzung von eva: eva EINGABEDATEI AUSGABEDATEI" << endl;
 cout << endl << "Aufbau der Eingabedatei: " << endl;
 cout << "nx ny       :  Anzahl Gitterpunkte in x- bzw. y- Richtung" << endl;
 cout << "Alpha       :  Einfallswinkel" << endl;
 cout << "Polarisation : =0 parallele Pol.  =1 senkrechte Pol." << endl;
 cout << "Anzahl Reflexionen" << endl;
 cout << "r0  l0      :  r0=Radius der Kugel  l0=Vakuumwellenlaenge" << endl;
 cout << "N0          :  Brechungsindex der Kugel" << endl;    
}
