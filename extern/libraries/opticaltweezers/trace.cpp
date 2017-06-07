#include <trace.h>
#include <surfArray.h>
extern ofstream oslog;  
#define CLIGHT 299792458 // m/s
#define MAX_RECURSIONS 20
#define DEBUG_INTENS
// #define DEBUG_RAY
int anzRaysFailed;
int anzRays;
//int hitStatus ;


int optFlag;
// surfArray SAF,SAL;

complex<double> ns=1.0; /// Brechungsindex Umgebungsmedium (Wasser)

int nP=0;
bool logP=false;
PointInfo *PI;

void saveForces(char *fn, int anzF, ForceData *F)
{
	ofstream os;
	os.open(fn);
	for (int i = 0; i < anzF; i++)
		os << F[i].P << "   " << F[i].F << endl;
	os.close();
}

void traceoneRay(Ray_pow &ray, Vector<double> *F, Vector<double> *L, complex<double> ns,int &iRC, int AnzReflex)
{
	/**
	   Berechnet die Kräfte auf die Partikel für einen einzigen Strahl (ray) 
	   F[Einschlussnr]=Kraft, die auf Einschluss <Einschlussnr> wirkt
	   P[Einschlussnr]=Leistung, die auf Einschluss <Einschlussnr> fällt
	   ns: Brechungsindex Umgebung
	*/

	Vector<double> Pe,P1,P2,er,r;
	Ray_pow tray;
	Vector<double> fe,fr,ft,fg;
	int iR;
	bool imEinschluss=false;
	bool Abbruch=false;
	complex<double> n1,n2;
	Vector<double> norm;
	iR=0;
	Pe=ray.P;
	anzRays++;
//	cout << "ray:" << ray << endl;
	if (iRC<MAX_RECURSIONS)
	{
	do
	{		
		imEinschluss=ray.imEinschluss;
		P1=ray.P;

		ray.next();          /// Nächsten Schnittpunkt suchen und Sprung durchführen		
		P2=ray.P;

	//	if (ray.imEinschluss || ray.einindex>-1)  cout << P1 << "  " << P2 << endl;
#ifdef DEBUG_RAY
		oslog << P1 << "     " << P2 << " % k=" << ray.k << endl;
#endif
		if (!ray.imEinschluss)  /// Strahl befindet sich nicht im Partikel
		{
		if (ray.einindex>-1)   /// Strahl trifft auf Partikel
		{
         n1=ns;              
		 n2=ray.Ein[ray.einindex]->getn();
         norm=ray.Ein[ray.einindex]->norm(ray.P);	  	// Oberflächennormale, muss entgegen dem einfallenden Strahl zeigen
		 fe=(ray.k*ray.Pow);//*abs(ray.k*norm);         // Kraftkomponente einfallender Strahl                
		 
		 tray=ray.reflect(norm,n1,n2);		           // Reflektiere Strahl 
		 fr=(-ray.k*ray.Pow);                          // Kraftkomponente reflektierter Strahl
		 ft=(-tray.k*tray.Pow)*real(n1/n2);//*abs(tray.k*norm); // Kraftkomponente gebrochener Strahl 
	     fg=fe+fr+ft;
	
		 if (optFlag & SAVE_FORCE)
		 {
			 if (tray.Ein[tray.einindex]->type == SURFACE)
			 {
	//			 SAF.addValue(tray.einindex, (surface *)(tray.Ein[tray.einindex]), ((surface *)(tray.Ein[tray.einindex]))->getCurrentIndex(), fg);
			 }
		 }

		 

		 r=ray.P-ray.Ein[tray.einindex]->P;            // Vektor vom Mittelpunkt des Objekts zum momentanen Punkt  
		 if (optFlag & SAVE_TORQUE)
		 {
			 if (tray.Ein[tray.einindex]->type == SURFACE)
			 {
//				 SAL.addValue(tray.einindex, (surface *)(tray.Ein[tray.einindex]), ((surface *)(tray.Ein[tray.einindex]))->getCurrentIndex(), r%fg);
			 }
		 }
		 F[tray.einindex]+=fg;
		 L[tray.einindex]+=r%fg;
		 if (optFlag && SAVE_FORCE) 
		 if (logP) 
		  {
							  if (nP==0) PI=(PointInfo *) malloc (sizeof(PointInfo) );
			  else PI=(PointInfo *) realloc (PI,sizeof(PointInfo) * (nP+1));
			  PI[nP].P=ray.P;
			  PI[nP].F=fg;
			  PI[nP].L=r%fg;
			  PI[nP].Pold=P1;
			  PI[nP].k=ray.k;
			  nP++;
		 } // if logP
		  ray=tray;	
		}
		else Abbruch=true;               /// Strahl hat kein Partikel getroffen
		}
		else   // Strahl befindet sich im Partikel
		{
			n2=ns;
		    n1=ray.Ein[ray.einindex]->getn();
			norm=-ray.Ein[ray.einindex]->norm(ray.P);            // Oberflächennormale zeigt nach Innen			
		    fe=(ray.k*ray.Pow)*real(n2/n1);//*abs(ray.k*norm);   // force because of incoming ray
		    tray=ray.reflect(norm,n1,n2);                        // tray->Strahl der herausgeht
		    fr=(-ray.k*ray.Pow);//*abs(ray.k*norm);              // Reflektierter Kraftanteil   
		    ft=(-tray.k*tray.Pow);//*abs(tray.k*norm);           // Transmittierter Kraftanteil 
#ifdef DEBUG_RAY
		oslog << "%% XXXkraus=" << tray.k << "     " << tray.P  << endl;
#endif
		fg=ft+fe+fr;
		if (optFlag & SAVE_FORCE)
		{
			if (ray.Ein[ray.einindex]->type == SURFACE)
			{
//				cout << "currentIndex=" << ((surface *)(ray.Ein[ray.einindex]))->getCurrentIndex() << " index=" << ray.einindex << endl;
	//			SAF.addValue(ray.einindex, (surface *)(ray.Ein[ray.einindex]), ((surface *)(ray.Ein[ray.einindex]))->getCurrentIndex(), fg);				
			}
		}


		F[ray.einindex]+=fg;
		r=ray.P-ray.Ein[ray.einindex]->P;  // Vektor vom Mittelpunkt des Objekts zum momentanen Punkt

		if (optFlag & SAVE_TORQUE)
		{
			if (ray.Ein[ray.einindex]->type == SURFACE)
			{
	//			SAL.addValue(ray.einindex, (surface *)(ray.Ein[ray.einindex]), ((surface *)(ray.Ein[ray.einindex]))->getCurrentIndex(), r%fg);
			}
		}


		L[ray.einindex]+=r%fg;
		if (logP)          // Nur bei Einzelrechnung: Damit an jeder Stelle die Strahlinfos abgespeichert werden können
		  {
			  if (nP==0) PI=(PointInfo *) malloc (sizeof(PointInfo) );
			  else PI=(PointInfo *) realloc (PI,sizeof(PointInfo) * (nP+1));
			  PI[nP].P=ray.P;
			  PI[nP].F=fg;
			  PI[nP].L=r%fg;
			  PI[nP].Pold=P1;
			  PI[nP].k=tray.k;
			  nP++;
		 } // if logP
		  iRC++;
		  traceoneRay(tray,F,L,ns,iRC,AnzReflex); //Rausgehender Strahl wird weiterverfolgt 
		  iR++;		  
		}		
	}
	while ((iR<AnzReflex) && (!Abbruch) && (ray.Pow>1E-20) && (abs(P1-P2)>1E-20));
	}
	else anzRaysFailed++;
    // cout << ray.P << " % Max. Anz. Rekursionen überschritten" << endl;
}

void trace (int nLS, LightSrc **ls, Vector<double> *F, Vector<double> *L)
{
	double alpha;
	double I;
	int statusLS;
	Ray_pow ray;
	Vector<double> P,k,P1,P2;
	Vector<double> **f = new Vector<double>*[nLS];  // Kräfte ohne c/n und ohne Rücksicht auf geeignete Normierung der Lichtquelle
	Vector<double> **l = new Vector<double>*[nLS];  // Drehmomente " "
	for (int i = 0; i < nLS; i++)
	{
		f[i] = new Vector<double>[ls[0]->AnzObjs];
		l[i] = new Vector<double>[ls[0]->AnzObjs];
	}
	//F[0]=zero;
	int c=0;
	
	anzRaysFailed=0;
	anzRays=0;
	for (int i=0; i<nLS; i++)
	{	 
	 c=0;
	 I = 0;
	 do
	 {
	    c++;	 
		statusLS=ls[i]->next(ray);       /// Nächsten Strahl holen		
		// cout << "wvl" << ls[i]->wvl << endl;
#ifdef DEBUG_RAY
		oslog << "%----------------------------------------------------------------" << endl;
		oslog << "% theta=" << asin(ray.P[0]/5E-6)/M_PI*180.0 << endl;
#endif

	//	cout << ray.P << "   " << ray.k << "  " << ray.r0 << endl ;
 		P=ray.P;
		k=ray.k;
		I=I+abs2(ray.E2);
		int iRC=0;
		traceoneRay(ray,f[i],l[i],ls[i]->n0,iRC);             /// Strahl verfolgen
		if (statusLS!=LIGHTSRC_ERROR) anzRaysFailed++;
	//	cout << P << "   " << k << "    " << F[0] << endl;
	//	F[0]=zero;
     } 
	 while ((statusLS!=LIGHTSRC_IS_LAST_RAY) /*&& (statusLS!=LIGHTSRC_ERROR)*/);  /// Schleife über alle Strahlen	 
	 for (int j = 0; j < ray.AnzEin; j++)
	 {
		 f[i][j] = f[i][j] / I*ls[i]->P0 * real(ls[0]->n0) / CLIGHT;
		 l[i][j] = l[i][j] * 1E-6 / I*ls[i]->P0 * real(ls[0]->n0) / CLIGHT; // 1E-6 weil alle Abstände in µm angegeben werden. l ist wird dann in Nm angegeben !
	 }

	}	
	
#ifdef MIT_NORMIERUNG	
	for (int j=0; j<ray.AnzEin; j++)  if (fabs(I)<1E-30) { F[j]=F[j]/I; L[j]=L[j]/I; }	
#else 
	for (int i = 0; i < nLS; i++)
	for (int j = 0; j < ray.AnzEin; j++) 
	{
		F[j] += f[i][j];
		L[j] += l[i][j];
	}
#endif
//	cout << "%" << anzRays << "   " << anzRaysFailed << "  %ANZRAYS" << endl;
	for (int i = 0; i < nLS; i++)
	{
		delete[] f[i];
		delete[] l[i];
	}
}
 

void findRootF(Vector<double> &r,int objIndex, int nLS, LightSrc **ls, Vector<double> *F,Vector<double> *L)
{
  // 
  for (int l=0; l<nLS; l++)
  {
   ls[l]->Ein[objIndex]->P=r;

  }
  
}

void doDynamics(int nLS, LightSrc **ls, int nSteps, double dt, double *rho, DynCalc **&Erg)
{
	Erg = new DynCalc*[ls[0]->AnzObjs];
	for (int i = 0; i < ls[0]->AnzObjs; i++)
		Erg[i] = new DynCalc[nSteps];

	Vector<double> *Q=new Vector<double> [ls[0]->AnzObjs];
    Vector<double> *L=new Vector<double> [ls[0]->AnzObjs];
    Matrix<double> *I=new Matrix<double> [ls[0]->AnzObjs];
	Matrix<double> *Iinv=new Matrix<double> [ls[0]->AnzObjs];
    Vector<double> w;
    Vector<double> n;
	Vector<double> *a=new Vector<double> [ls[0]->AnzObjs];
	Vector<double> *v=new Vector<double> [ls[0]->AnzObjs];
    double phi,absw;
    Matrix<double> *D=new Matrix<double> [ls[0]->AnzObjs];
	double *m=new double[ls[0]->AnzObjs];
	int anzObjs=ls[0]->AnzObjs;
	// Erst mal Trägheitstensoren berechnen 
	// Trägheitsmoment wird 
	// double rho=1.18E-15; // Dichte (PMMA) in kg/µm³  
    for (int j=0; j<anzObjs; j++)
	{		
		D[j]=ls[0]->Ein[j]->H;               // Orientierung übernehmen
		m[j]=ls[0]->Ein[j]->Volume()*rho[j]; // Achtung Volumen wird in µm³ angeben !!!
		I[j]=computeInertia(ls[0]->getObject(j) ) * m[j]; // Trägheitstensor berechnen in kg*m²
		Iinv[j]=invert(I[j]);               // Es wird eigentlich die Inverse vom Trägheitstensor gebraucht, denn 
		                                    // L=I*w --> w=I^-1*L
	}

	for (int i=0; i<nSteps; i++)
	{
        // Erst mal Kräfte und Drehmomente berechnen 
		for (int l=0; l<nLS; l++) { Q[l]=zero; L[l]=zero; ls[l]->reset(); }
        trace(nLS,ls,Q,L);

		
        for (int j=0; j<anzObjs; j++)
		{
			// Neuen Ort berechnen
			// Kraft Q=F/P*c/n => F=Q
			a[j] = Q[j]  / m[j]; // Beschleunigung in m/s^2 
			v[j] += a[j] * dt;   // Geschwindigkeit in m/s                                      // Geschwindigkeit in µm/s;					
			w = Iinv[j] * L[j] * dt; // Winkelgeschwindigkeit in 1/s
				
			absw=abs(w);
			n=w/absw;    // Rotationsachse
			phi=absw*dt; // Änderung des Drehwinkels
	//		cout << i << "   " << ls[0]->Ein[j]->P << "   " << n << "   " << phi/M_PI*180.0 << endl;
			Erg[j][i].angle = phi;
			Erg[j][i].axis = n;
			Erg[j][i].Pos = ls[0]->Ein[j]->P;

			D[j]*=drehmatrix(n,phi); // 
			for (int l=0; l<nLS; l++) 
			{ 
				ls[l]->Ein[j]->P+=v[j]*dt*1E6;          // Ort in µm
				ls[l]->Ein[j]->setMatrix(D[j]);     // Orientierung im Raum ändern
			}

			// cout << ls[0]->Ein[j]->P << "  " << v[j] << "  " << Q[j] << "  " << ls[0]->Ein[j]->e[0] << "  " << ls[0]->Ein[j]->e[1] << "  " << ls[0]->Ein[j]->e[2] << endl; 
		//	cout << ls[0]->Ein[j]->P <<  "  " << ls[0]->Ein[j]->e[0] << "  " << ls[0]->Ein[j]->e[1] << "  " << ls[0]->Ein[j]->e[2] << "   " << w << endl; 
  
		}

	}
	delete[] Q;
	delete[] L;
	delete[] I;
	delete[] Iinv;
	delete[] a;
	delete[] v;
	delete[] D;
	delete[] m;
}

int iRC;
void traceoneRay(Ray_pow &ray, int &anzP, Vector<double> **&P, int &anzF, ForceData *&fd, complex<double> ns, int reflexCounter, int AnzReflex, bool drawOutgoing)
{
	/**
	Berechnet die Schnittpunkte der Strahlen bzw. Strahlabschnitte mit den Partikeloberflächen
	P[i,j]: Start- (i=0) und Endpunkte (i=1) der Strahlabschnitte
	ns: Brechungsindex Umgebung
	*/
	bool hatGetroffen=false;
	Vector<double> P1,P2,Pe,  er, r;
	Ray_pow tray,hray;
	int iR;
	bool imEinschluss = false;
	bool Abbruch = false;
	complex<double> n1, n2;
	Vector<double> norm;
	Vector<double> fr, fe, fg, ft;
	bool saveRay;
	iR = 0;
	Pe = ray.P;
	anzRays++;
	iRC++;

	if (iRC < MAX_RECURSIONS)
	{
		do
		{
			imEinschluss = ray.imEinschluss;
			
			
			P1 = ray.P;
			ray.next();          /// Nächsten Schnittpunkt suchen und Sprung durchführen	
			P2 = ray.P;
//			cout << P1 << "   " << P2 << endl;
		
#ifdef DEBUG_RAY
			oslog << P1 << "     " << P2 << " % k=" << ray.k << endl;
#endif
			if (!ray.imEinschluss)  /// Strahl befindet sich nicht im Partikel
			{
				if (ray.einindex > -1)   /// Strahl trifft auf Partikel
				{
					P = (Vector<double> **) realloc(P, sizeof(Vector<double> *) * (anzP + 1));
					P[anzP] = (Vector<double> *) malloc(sizeof(Vector<double>) * 2);
					P[anzP][0] = P1;
					P[anzP][1] = P2;
					anzP++;

					n1 = ns;
					n2 = ray.Ein[ray.einindex]->getn();
					norm = -ray.Ein[ray.einindex]->norm(ray.P);
					fe = (ray.k*ray.Pow);					         // Kraftkomponente einfallender Strahl                

					r = ray.P - ray.Ein[tray.einindex]->P;
					tray = ray.reflect(-norm, n1, n2);
					
					if (saveRay)
					{
						fr = (-ray.k*ray.Pow);                           // Kraftkomponente reflektierter Strahl
						ft = (-tray.k*tray.Pow)*real(n1 / n2);           // Kraftkomponente gebrochener Strahl 
						fg = fe + fr + ft;
					}
					
					traceoneRay(tray, anzP, P, anzF, fd, ns, reflexCounter, AnzReflex,drawOutgoing); //Rausgehender Strahl wird weiterverfolgt 					
					reflexCounter++;
				}
				else
				{
				Abbruch = true;               /// Strahl hat kein Partikel getroffen
				if (drawOutgoing)
				{
					P = (Vector<double> **) realloc(P, sizeof(Vector<double> *) * (anzP + 1));
					P[anzP] = (Vector<double> *) malloc(sizeof(Vector<double>) * 2);
					P[anzP][0] = P1;
					P[anzP][1] = P2;
					anzP++;
				}
				}
			}
			else   // Strahl befindet sich im Partikel
			{
				P = (Vector<double> **) realloc(P, sizeof(Vector<double> *) * (anzP + 1));
				P[anzP] = (Vector<double> *) malloc(sizeof(Vector<double>) * 2);
				P[anzP][0] = P1;
				P[anzP][1] = P2;
				anzP++;
				hatGetroffen = true; 
				n2 = ns;
				n1 = ray.Ein[ray.einindex]->getn();
				norm = -ray.Ein[ray.einindex]->norm(ray.P);              // Oberflächennormale zeigt nach Innen			
				fe = (ray.k*ray.Pow)*real(n2 / n1);					     // force because of incoming ray
				tray = ray.reflect(norm, n1, n2);                        // tray->Strahl der herausgeht
				if (saveRay)
				{
					fr = (-ray.k*ray.Pow);                                   // Reflektierter Kraftanteil   
					ft = (-tray.k*tray.Pow)*real(n1 / n2);		             // Transmittierter Kraftanteil 
					fg = ft + fe + fr;
					anzF++;
				}
				traceoneRay(tray, anzP, P, anzF, fd, ns, reflexCounter, AnzReflex,drawOutgoing);      // Rausgehender Strahl wird weiterverfolgt 				
				reflexCounter++;
			}
		} while ((reflexCounter<AnzReflex)  && (!Abbruch) );
	}
	//else anzRaysFailed++;
	
	// cout << ray.P << " % Max. Anz. Rekursionen überschritten" << endl;
}
 void trace(int nLS, LightSrc **ls, int &anzP, Vector<double> **&P, int &anzF, ForceData *&F, int AnzReflex , bool drawOutgoing)
/**
Strahlverfolgungsroutine für Strahldarstellung (Berechnet nur die Haltepositionen/Schnittpunkte mit Oberflächen
nLS : Anzahl Strahlquellen
ls : Strahlquellen
Rückgabeparameter:
anzP : Anzahl Punkte 
P : Punkte 
*/
{
	double alpha;
	double I;
	int statusLS;
//	int iRC;
	Ray_pow ray;
	Vector<double>  k;
	P = (Vector<double> **) malloc(sizeof(Vector<double> *));
	anzP = 0;
	anzF = 0;
	int c = 0;
#ifdef DEBUG_INTENS
	ofstream os;
	os.open("intensity.dat");
#endif
	cout << "D=" << ls[0]->D << "    dens=" << ls[0]->density << "   Anz=" << ls[0]->N << endl;
	anzRaysFailed = 0;
	anzRays = 0;
	for (int i = 0; i < nLS; i++)
	{
		c = 0;
		do
		{
			c++;
			statusLS = ls[i]->next(ray);       /// Nächsten Strahl holen				
#ifdef DEBUG_INTENS
			os << ray.E1 << "   " << ray.E2 << "  " << ray.Pow << endl;
#endif
			iRC = 0;
			traceoneRay(ray, anzP, P, anzF, F, ls[i]->n0, 0,AnzReflex,drawOutgoing);             /// Strahl verfolgen
			if (statusLS != LIGHTSRC_ERROR) anzRaysFailed++;
			
		} while (statusLS != LIGHTSRC_IS_LAST_RAY);
	}
//	cout << "%" << anzRays << "   " << anzRaysFailed << "  %ANZRAYS" << endl;
	// saveForces("test1.dat", anzP, F); 
#ifdef DEBUG_INTENS
	os.close();
#endif
}
