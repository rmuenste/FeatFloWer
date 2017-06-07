#include "LightSrc.h"
#include "PStrahl.h"
#include "misc.h"

double erf(double x)
{    
	// constants    
	double a1 =  0.254829592;    
	double a2 = -0.284496736;    
	double a3 =  1.421413741;    
	double a4 = -1.453152027;    
	double a5 =  1.061405429;    
	double p  =  0.3275911;    
	
	// Save the sign of x    
	int sign = 1;    
	if (x < 0)        
		sign = -1;    
	x = fabs(x);    // A&S formula 7.1.26    
	double t = 1.0/(1.0 + p*x);    
	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);    
return sign*y;
}

LightSrc::LightSrc(const LightSrc &L)
{
 Pos=L.Pos;
 type=L.type;
 density=L.density;
 e1=L.e1;
 e2=L.e2;
 k=L.k;
 i1=L.i1;
 i2=L.i2;
 Pol=L.Pol;
 r0=L.r0;
 wvl=L.wvl;
 P0 = L.P0;
 AnzObjs=L.AnzObjs;
// Ein=L.Ein;
 copyFormList(Ein,L.Ein,AnzObjs);
 D=L.D;
 raytype=L.raytype; 
 n0=1.0;
 N=L.N;
}

void LightSrc::binRead (ifstream &is)
{
	double nre,nim;
	Pos.binRead(is);
	is.read((char *) &type,(char) sizeof(type));
	is.read((char *) &density,(char) sizeof(density));
	k.binRead(is);
	is.read((char *) &N,(char) sizeof(N));
	Pol.binRead(is);
	is.read((char *) &r0,(char) sizeof(r0));
	is.read((char *) &wvl,(char) sizeof(wvl));
	is.read((char *) &AnzObjs,(char) sizeof(AnzObjs));
	is.read((char *) &nre, (char) sizeof(nre));
	is.read((char *) &nim, (char) sizeof(nim));
	n0=complex<double> (nre,nim);
	is.read((char *) &D, (char) sizeof(D));
	e1.binRead(is);
	e2.binRead(is);
	is.read((char *) &raytype, (char) sizeof(raytype));
	switch (type)
	{
	 case LIGHTSRC_SRCTYPE_PLANE : ((LightSrcPlane *) this)->binReadItem(is); break;
	 case LIGHTSRC_SRCTYPE_GAUSS : ((LightSrcGauss *) this)->binReadItem(is); break;
	}
}

void binWriteLSList (ofstream &os, int nLS, LightSrc **ls)
{
	for (int i=0; i<nLS; i++)
	{
		os.write((char *) &ls[i]->type,sizeof(ls[i]->type));
		ls[i]->binWrite(os);
		switch	(ls[i]->type)
	    {
	     case LIGHTSRC_SRCTYPE_PLANE : ((LightSrcPlane *) ls[i])->binWriteItem(os); break;
	     case LIGHTSRC_SRCTYPE_GAUSS : ((LightSrcGauss *) ls[i])->binWriteItem(os); break;
	    }
	}
}

void binReadLSList (ifstream &is, int nLS, LightSrc ** &ls)
{
	int type;
	ls=(LightSrc **) malloc (sizeof (LightSrc *) * nLS);
	for (int i=0; i<nLS; i++)
	{
		
		is.read ((char *) &type, sizeof(type));
		switch	(type)
	    {
	     case LIGHTSRC_SRCTYPE_PLANE : 
			 ls[i]=new LightSrcPlane();
			 ls[i]->binRead(is);
			 ((LightSrcPlane *) ls[i])->binReadItem(is); 
			 break;
	     case LIGHTSRC_SRCTYPE_GAUSS : 
			 ls[i]=new LightSrcGauss();
			 ((LightSrcGauss *) ls[i])->binRead(is);
			 ((LightSrcGauss *) ls[i])->binReadItem(is); 
			 break;
	    }
	}
}


void LightSrc::binWrite (ofstream &os)
{
	double nre,nim;
	Pos.binWrite(os);
	os.write((char *) &type,(char) sizeof(type));
	os.write((char *) &density,(char) sizeof(density));
	k.binWrite(os);
	os.write((char *) &N,(char) sizeof(N));
	Pol.binWrite(os);
	os.write((char *) &r0,(char) sizeof(r0));
	os.write((char *) &wvl,(char) sizeof(wvl));
	os.write((char *) &AnzObjs,(char) sizeof(AnzObjs));
	nre=real(n0);
	nim=imag(n0);
	os.write((char *) &nre, (char) sizeof(nre));
	os.write((char *) &nim, (char) sizeof(nim));
	os.write((char *) &D, (char) sizeof(D));
	e1.binWrite(os);
	e2.binWrite(os);
	os.write((char *) &raytype, (char) sizeof(raytype));	
	switch (type)
	{
	 case LIGHTSRC_SRCTYPE_PLANE : ((LightSrcPlane *) this)->binWriteItem(os); break;
	 case LIGHTSRC_SRCTYPE_GAUSS : ((LightSrcGauss *) this)->binWriteItem(os); break;
	}

}



void LightSrc::ObjectList (int Anz, Form **Ein) 
{ 
	AnzObjs=Anz; 
	copyFormList(this->Ein,Ein,Anz);
//	this->Ein=Ein; 
}

void LightSrc::setPos(Vector<double> P) 
{ 

	Pos = P;
}

void LightSrc::setObject(Form *O, int i)
/**
   O : Zeiger auf das Objekt
   i : Index des Objekts, das geändert werden soll
       ist i<0 oder größer als die Anzahl von vorhandenen Objekten, dann wird das Objekt als
	   neues Element an das Ende angehängt
**/
{
  if ((i<0) || (i>AnzObjs-1)) 
  {
	Ein=(Form **) realloc (Ein,sizeof (Form *) * (AnzObjs+1));
	copyInc(Ein[AnzObjs],O);
	AnzObjs++;
  }
  else
  {	  
   delete Ein[i];
   	copyInc(Ein[i],O);
  }
}

LightSrcPlane::LightSrcPlane(void)
{
	type=LIGHTSRC_SRCTYPE_PLANE;
	k=ez;
	density=1;
	raytype=LIGHTSRC_RAYTYPE_IRAY;
	AnzObjs=0;
	Ein=0;
	e1=ex;
	e2=ey;

	reset();
}



LightSrcPlane::LightSrcPlane (Vector<double> Pos, int N, double wvl, double r0, double D, Vector<complex<double> > Pol, int raytype)
{
	e1=ex;
	e2=ey;
	
 this->Pos=Pos;
 this->density=D/((double) N);
 this->type=LIGHTSRC_SRCTYPE_PLANE;
 this->k=ez;
 this->raytype=raytype;
 this->Pol=Pol;
 this->r0=r0;
 this->wvl=wvl;
 this->D=D; 
 this->N=N;
 AnzObjs=0;
 Ein=0;
 reset();
}

LightSrcPlane::LightSrcPlane (const LightSrcPlane &)
{
  type=LIGHTSRC_SRCTYPE_PLANE;
}

int LightSrcPlane::next(IStrahl &S)
{  	 
	
 Ebene E;
  Vector<double> P=Pos+(i1*density)*e1+(i2*density)*e2;
  E.e1=e1;
  E.e2=e2;
  E.n=k;
//	Vector<double> P(0.4,0.0,-2.0);
  S=IStrahl(P,Pol,k,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);  
  S.E1=Pol;
  S.E2=Pol;
  S.init_Efeld(E,Pol);
  i1++;
   
  if (i1*density>D) {i1=0; i2++; cout << "i2=" << i2 << endl; }
  if (i2*density>=D) {return LIGHTSRC_IS_LAST_RAY; }
  return LIGHTSRC_NOT_LAST_RAY;
}

int LightSrcPlane::next(Ray_pow &S)
{  	 
	
 Ebene E;
 double Pow;
 
 Vector<double> P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2; 
// Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
  E.e1=e1;
  E.e2=e2;
  E.n=k;
  Pow=1.0/((double)(N*N)*D*D);
  S=Ray_pow(Pow,P,Pol,k,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);  
  S.init_Efeld(E,Pol);
  S.P=P;
  S.E1=Pol;
  S.E2=sqrt(Pow)*Pol;
  S.k=k;  
  i1++;
   
 // if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
  if (i1*density>D) {i1=0; i2++; }
  if (i2*density>=D) {return LIGHTSRC_IS_LAST_RAY; }
  return LIGHTSRC_NOT_LAST_RAY;
}


int LightSrcPlane::next(Ray &S)
{  	 	
   Vector<double> P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2;  
//	Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
//	cout << "NEXT: P=" << P << endl;
  S=Ray(P,density,density,Pol,k,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);  
  S.setN0(n0);
  // S.init_Efeld(Pol,1);
  i1++;
//  if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!
  if (i1*density>D) {i1=0; i2++; }
  if (i2*density>=D) return LIGHTSRC_IS_LAST_RAY;
  return LIGHTSRC_NOT_LAST_RAY;
}

void LightSrc::reset()
{
 i1=0;
 i2=0; 
}


LightSrc::~LightSrc(void)
{
 if (AnzObjs>0) delete[] Ein;
 AnzObjs=0;
 reset();
}

void LightSrc::clearObjects()
{
 if (AnzObjs>0) 
 {
	 for (int i=0; i<AnzObjs; i++)
		 delete Ein[i];
	 delete[] Ein;
	 AnzObjs=0;
 }
}

LightSrcGauss::LightSrcGauss(void)
{
	type=LIGHTSRC_SRCTYPE_GAUSS;
	k=ez;
	density=1;
	raytype=LIGHTSRC_RAYTYPE_IRAY;
	AnzObjs=0;
	Ein=0;
	e1=ex;
	e2=ey;
	n0=1.0;
	Pol=Vector<complex<double> > (0,1.0,0);
	polType=LIGHTSRC_POL_Y;
	reset();
}


LightSrcGauss:: LightSrcGauss (Vector<double> Pos,  int N, double wvl, double w0, Vector<double> focuspos, double D,  Vector<complex<double> > Pol, int raytype, double r0)
/**
  Konstruktor für Gauss-Strahlungsquelle 
  Strahlen laufen aus einem Rechteck am Ort Pos der Breite D in z-Richtung auf den Fokuspunkt zu
  Parameter:
  Pos : Position der Quelle (Mitte) 
  N   : Anzahl Strahlen je Raumrichtung
  wvl : Wellenlänge
  w0  : Taillendurchmesser
  focuspos : Position des Fokus
  D   : Breite der Lichtquelle 
  r0  : Radius der "Weltkugel"
**/
{
	/*
  	e1=ex;
	e2=ey;
	*/
	// Erst mal die Vektoren e1 und e2 setzen
	
 this->Pos=Pos;
 this->density=D/((double) N);
 this->type=LIGHTSRC_SRCTYPE_GAUSS;
 k=focuspos-Pos;
 k=k/abs(k);
 e1=k%ez;
 if (abs(e1)<1E-10) e1=ex;
 e2=k%e1;
 e1=e1/abs(e1);
 e2=e2/abs(e2);

 this->raytype=raytype;
 this->Pol=Pol;
 this->r0=r0;
 this->wvl=wvl;
 this->D=D;  
 this->focuspos=focuspos;
 Pol=Vector<complex<double> > (0,1.0,0);
	polType=LIGHTSRC_POL_Y;

 this->Pol=Pol; 
 this->w0=w0;
 this->N=N;
 this->P0=P0;
 f=abs(Pos-focuspos);
 z0=M_PI*w0*w0/wvl;
 AnzObjs=0;
 Ein=0; 
 n0=1.0;
 double d=abs(Pos-focuspos);
 calcz0();
 calcw(d);
 double theta = atan(wvl / (M_PI*w0));
 NA = real(n0)*sin(theta);
 reset();
}


void LightSrcGauss::setWvl(double wvl)
{
	this->wvl = wvl;
	double theta = asin(real(NA/this->n0)); 
	w0 = wvl / (M_PI*tan(theta));
}

int LightSrcGauss::next(IStrahl &S)
{  	 
    Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2; 
	 // Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!! 
 double x1,x2,x3;
 x1=P*e1;
 x2=P*e2;
 double r2=x1*x1+x2*x2;
  
 double s2=w0*w0/log(2.0);
 double g;
 double L=0.1;
 double R;
 
 Vector<double> h,hk; 
 Matrix<double> DM;
 Vector<complex<double> > E;
 complex<double> E0;
 double absh,gamma;
  double zeta;
  double absfp;

  fp=focuspos-P;  // Hier beginnt der Strahl
  hk=fp/abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus
  
 
  h=k%hk;          
  absh=abs(h);

  x3=fp*k;
  R=x3*(1.0+z0*z0/(x3*x3));
  if (z0==0) zeta=M_PI/2.0;
  else zeta=atan(x3/z0);

  
  if (absh==0) 
  {
	  hk=k; 
	  gamma=0;
  }
  else 
  {
  h/=absh;  
  gamma=acos(k*hk/(abs(k)*abs(hk)));    
  }
 
  Pol=Vector<complex<double> >(0,1,0);  // y-Polarisation
  
  DM=drehmatrix(h,gamma); 
  S=IStrahl(P,Pol,hk,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);
  S.k=k;  
  S.n=n0;
  absfp=fp*k;
   g=gaussw(-fp*k,2.0*M_PI/real(S.k0),w0); // w(z)   
r2=S.P[0]*S.P[0]+S.P[1]*S.P[1];
 
//E0=sqrt(1.0/M_PI*g/(absfp*w0*w0)*sqrt(2.0/M_PI)*exp(-2.0*absfp*absfp/(g*g))
//		*1.0/(1.0-erf(sqrt(2.0)*absfp/g)))*w0/g*exp(-r2/g/g);

  E0=sqrt(P0)*sqrt(2.0/M_PI/g/g)*exp(-r2/g/g); // ohne Korrektur
// E0=1;
Vector<double> F=focuspos;

     S.k=focuspos-S.P;   // Richtungsvektor auf den Fokus gerichtet
	 S.k/=abs(S.k);	      
	 h=k%S.k;               
	 absh=abs(h);	
	 if (absh!=0) 
	 {
		 h/=absh;
	     gamma=acos(S.k*k);     
         DM=drehmatrix(h,gamma);
         S.E1=E0*DM*Pol;	 	
	 }
	 else S.E1=E0*Pol;
	 S.n=n0;
  i1++;  
  if (i1*density>=D) {i1=0; i2++; }
  if (i2*density>=D) return LIGHTSRC_IS_LAST_RAY;
  return LIGHTSRC_NOT_LAST_RAY;   
}

int LightSrcGauss::next(Ray_pow &S)
{  	 
    Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2; 
//  Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1; // nur zu TESTZWECKEN !!!
	// P : Startort, density=Anzahl Strahlen/Längeneinheit, D: Breite der Lichtquelle
	//     i1=momentaner Index in e1-Richtung (normalerweise x-Richtung), i2=momentaner Index in e2-Richtung (normalerweise y-Richtung)
	
 double x1,x2,x3; // Hilfsgrößen
 x1=P*e1;
 x2=P*e2;
 double r2=x1*x1+x2*x2;  // Quadrat des Abstands von der Laserstrahlachse
  
 double s2=w0*w0/log(2.0);  // w0 : Strahltaille 
 double g;            
 double L=0.1;
 double R;
 
 Vector<double> h,hk; 
 Matrix<double> DM;
 Vector<complex<double> > E;
 complex<double> E0;
 double absh,gamma;
  double zeta;
  double absfp;
  calcNormfak();
  fp=focuspos-P;  // Hier beginnt der Strahl
  hk=fp/abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus
  
 
  h=k%hk;          
  absh=abs(h);

  x3=fp*k;
  R=x3*(1.0+z0*z0/(x3*x3));
  if (z0==0) zeta=M_PI/2.0;
  else zeta=atan(x3/z0);

  
  if (absh==0) 
  {
	  hk=k; 
	  gamma=0;
  }
  else 
  {
  h/=absh;  
  gamma=acos(k*hk/(abs(k)*abs(hk)));    
  }
 
  // Pol=Vector<complex<double> >(0,1,0);  // y-Polarisation
  
  DM=drehmatrix(h,gamma); 
  S=Ray_pow(1,P,Pol,hk,n0,r0,2.0*M_PI/wvl,AnzObjs,Ein);
  S.k=k;  
  S.n=n0;
  absfp=fp*k;
  g=gaussw(-fp*k,2.0*M_PI/real(S.k0),w0); // w(z)  

 
  r2=S.P[0]*S.P[0]+S.P[1]*S.P[1];

  S.Pow=P0*k*hk*exp(-2.0*r2/g/g)/((double) N*N)*D*D/g;//*sqrt(8.0/M_PI);// *real(Normfak);
 // cout << "%w=" << g << "   Pow=" << S.Pow << "   P0=" << P0 << endl;
#ifdef MIT_NORMIERUNG
  S.Pow = k*hk*exp(-2.0*r2 / g / g) / ((double)N*N)*D*D / g;
#endif
  S.E1=Pol;
  S.E2=sqrt(S.Pow)*Pol;
  
Vector<double> F=focuspos;
     S.k=focuspos-S.P;   // Richtungsvektor auf den Fokus gerichtet
	 S.k/=abs(S.k);	      
	 h=k%S.k;               
	 absh=abs(h);	
	 if (absh!=0) 
	 {
		 h/=absh;
	     gamma=acos(S.k*k);     
         DM=drehmatrix(h,gamma);
		 S.E1=DM*S.E1;	 	
	     S.E2=DM*S.E2;	 	
	 }
	 // else S.E1=E0*Pol;
	 S.n=n0;
  i1++;  
  if (type==LIGHTSRC_SRCTYPE_TOPHAT_FOKUS) 
  {
	  S.E2=S.E2/abs(S.E2);
      S.Pow=1.0;
  }
  
//  if (i1*density>D) { return LIGHTSRC_IS_LAST_RAY; } // nur zu TESTZWECKEN !!!!
  if (i1*density>D) {i1=0; i2++;} 
  if (i2*density>D) return LIGHTSRC_IS_LAST_RAY;
  return LIGHTSRC_NOT_LAST_RAY;   
}



int LightSrcGauss::next(Ray &S)
{
	Vector<double> Ph=(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2; 
	Vector<double> P=Pos+Ph;
	Vector<double> k=focuspos-P;  // Richtung des Strahles
	k=k/abs(k);
	double r=abs(Ph);
	double s=w0/(2.0*sqrt(log(2.0))); 
	double E0=exp(-r*r/s/s);
	S=Ray(P,density,density,Pol,ez,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);
    S.setN0(n0);
    S.n=n0;
	for (int i=0; i<5; i++) S.E[i]=Pol*E0;
	i1++;
	 if (i1*density>D) {i1=0; i2++; }
  if (i2*density>=D) return LIGHTSRC_IS_LAST_RAY;
  return LIGHTSRC_NOT_LAST_RAY;
}

/*int LightSrcGauss::next(ray &S)
{  	 
    Vector<double> fp,P=Pos+(i1*density-D/2.0)*e1+(i2*density-D/2.0)*e2; 
	 // Vector<double> P=Pos+(i1*density-D/2.0)*e1; // NUR ZU TESTZWECKEN !!!!!!!!!! 
 double x1,x2,x3;
 x1=P*e1;
 x2=P*e2;
 double r2=x1*x1+x2*x2;
  
 double s2=w0*w0/log(2.0);
 double g;
 double L=0.1;
 double R;
 
 Vector<double> h,hk; 
 Matrix<double> DM;
 Vector<complex<double> > E;
 complex<double> E0;
 double absh,gamma;
  double zeta;
  double Eph;
  double l;
  double absfp;

  fp=focuspos-P;  // Hier beginnt der Strahl
  hk=fp/abs(fp);  // Normierter Richtungsvektor vom Startpunkt zum Fokus
  
 
  h=k%hk;          
  absh=abs(h);

  x3=fp*k;
  R=x3*(1.0+z0*z0/(x3*x3));
  if (z0==0) zeta=M_PI/2.0;
  else zeta=atan(x3/z0);

  
  if (absh==0) 
  {
	  hk=k; 
	  gamma=0;
  }
  else 
  {
  h/=absh;  
  gamma=acos(k*hk/(abs(k)*abs(hk)));    
  }
 
  Pol=Vector<complex<double> >(0,1,0);  // y-Polarisation
  
  DM=drehmatrix(h,gamma); 
  S=ray(P,density,density,Pol,hk,1.0,r0,2.0*M_PI/wvl,AnzObjs,Ein);
  S.setN0(n0);
  S.n=n0;
  absfp=fp*k;
  g=w0/(2.0*sqrt(log(2.0)));
  
   g=gaussw(-fp*k,2.0*M_PI/S.k0,w0); // w(z)   
r2=S.P[4][0]*S.P[4][0]+S.P[4][1]*S.P[4][1];
 
E0=sqrt(1.0/M_PI*g/(absfp*w0*w0)*sqrt(2.0/M_PI)*exp(-2.0*absfp*absfp/(g*g))
		*1.0/(1.0-erf(sqrt(2.0)*absfp/g)))*w0/g*exp(-r2/g/g);
		
//  E0=sqrt(2.0/M_PI/g/g)*exp(-r2/g/g); // ohne Korrektur
// E0=1;
Vector<double> F=focuspos;
double B,C2;

for (int i=0; i<5; i++) 
{
     S.k[i]=focuspos-S.P[i];   // Richtungsvektor auf den Fokus gerichtet
	 S.k[i]/=abs(S.k[i]);	      
	 h=k%S.k[i];               
	 absh=abs(h);	
	 if (absh!=0) 
	 {
		 h/=absh;
	     gamma=acos(S.k[i]*k);     
         DM=drehmatrix(h,gamma);
         S.E[i]=E0*DM*Pol;	 	
	 }
	 else S.E[i]=E0*Pol;
	 S.n0=n0;
}
  i1++;  
//  cout << "E=" << S.E[4] << endl;
 // if (i1*density>D) return LIGHTSRC_IS_LAST_RAY; // NUR ZU TESTZWECKEN !!!!!!!!!!!!!!!!!!
  if (i1*density>D) {i1=0; i2++; }
  if (i2*density>=D) return LIGHTSRC_IS_LAST_RAY;
  return LIGHTSRC_NOT_LAST_RAY;
   
}*/

LightSrcGauss::LightSrcGauss(const LightSrcGauss &L) : LightSrc(L)
{
 focuspos=L.focuspos;
 n0=L.n0;
 wvl = L.wvl;
 w0 = L.w0;
 P0 = L.P0;
 k = L.k;
 AnzObjs = L.AnzObjs;
}

void LightSrcGauss::initGauss(Vector<double> &P, Vector<double> &n, Vector<complex<double> > &E)
{

  Pol=E;
  k=n;
  focuspos=P;
  calcz0();
   double d=abs(Pos-focuspos);
  calcw(d);
  /*Vector<double> h,d;
  Matrix<double> DM;
  double r,absh,w2,Z,gamma;
  // double l,R,zeta;
  
  // double k0=2.0*M_PI/wvl;  
  h=focuspos-P;
  absh=abs(h);
  k=h/absh;   
  r=abs(h%this->k);    */

 
 // E=Vector<complex<double> >(0,E0,0);
  
}  
 

/*LightSrcGauss::LightSrcGauss(const LightSrc &L) : LightSrc(L) 
{
 
};*/

ostream& operator << (ostream &os, LightSrc *ls)
{
	os << "% Anzahl Objekte:" << ls->AnzObjs << endl;
	for (int i=0; i<ls->AnzObjs; i++)
	{		
		os << "&------------------- Objekt " << i << " -----------------" << endl;
		os << "& P=" << ls->Ein[i]->P << endl;
		if (ls->Ein[i]->type==ELLIPSOID)
		{
			FormEllipsoid *e=(FormEllipsoid *)ls->Ein[i];
			os <<  "% Halbachsen= " << e->r << endl;			
		}
		os << "% Brechungsindex n=" << real(ls->Ein[i]->n) << "+ i*" << imag(ls->Ein[i]->n) << endl;
	}

	os << "% D=" << ls->D << "   k=" << ls->getk() << "  N=" << ls->N << endl;
	os << "% density=" << ls->density << endl;
	return os;
}


void LightSrcGauss::binWriteItem(ofstream &os)
{		  
	os.write((char *) &w0, (char) sizeof(w0));
	focuspos.binWrite(os);
	os.write((char *) &z0, (char) sizeof(z0));
	os.write((char *) &P0, (char) sizeof(P0));
	os.write((char *) &f, (char) sizeof(f));	
}

void LightSrcGauss::binReadItem(ifstream &is)
{
	is.read((char *) &w0, (char) sizeof(w0));
	focuspos.binRead(is);
	is.read((char *) &z0, (char) sizeof(z0));
	is.read((char *) &P0, (char) sizeof(P0));
	is.read((char *) &f, (char) sizeof(f));
}

void LightSrcGauss::setNA(double NA)
{
	this->NA = NA; 
	double theta = asin(NA);
	double l = abs(Pos - focuspos);
	w0 = wvl / (M_PI*tan(theta));
	D = 2.0*l*wvl / (M_PI*w0);
	density = D / ((double)N);
	calcz0();
}

void copyLightSrcList (LightSrc **&d, LightSrc **s, int anz)
{
	if ((anz>0) && (s!=0))
    {
      d=(LightSrc **) malloc (sizeof (LightSrc *) * anz);
      for (int i=0; i<anz; i++)
	  {
		 switch (s[i]->type)
		{
		  case LIGHTSRC_SRCTYPE_PLANE : d[i]=new LightSrcPlane (*((LightSrcPlane *)s[i])); break;
		  case LIGHTSRC_SRCTYPE_GAUSS : d[i]=new LightSrcGauss (*((LightSrcGauss *)s[i])); break;
          case LIGHTSRC_SRCTYPE_TOPHAT_FOKUS : d[i]=new LightSrcGauss (*((LightSrcGauss *)s[i])); 
			                                   d[i]->type=LIGHTSRC_SRCTYPE_TOPHAT_FOKUS;
												break;
		 }
	  }
    }
}
