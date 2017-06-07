#ifndef Lightsrc_H
#define Lightsrc_H
#define LIGHTSRC_RAYTYPE_RAY 1    ///< Strahltypklasse : Strahl
#define LIGHTSRC_RAYTYPE_IRAY 2   ///< Strahltypklasse : IStrahl
#define LIGHTSRC_RAYTYPE_PRAY 3   ///< Strahltypklasse : PStrahl
#define LIGHTSRC_SRCTYPE_PLANE 1  ///< Lichtquelle ist eine ebene Welle
#define LIGHTSRC_SRCTYPE_GAUSS 2  ///< Lichtquelle ist ein Gaussstrahl
#define LIGHTSRC_SRCTYPE_TOPHAT_FOKUS 3 ///< Lichtquelle ist ein fokussierter Tophat
#include "vector.h"
#include "strahl.h"
#include "istrahl.h"
#include "form.h"
#include "PStrahl.h"
#define LIGHTSRC_NOT_LAST_RAY 0  ///< Der gerade betrachtete Strahl ist nicht der letzte Strahl 
#define LIGHTSRC_IS_LAST_RAY 1   ///< Der gerade betrachtete Strahl ist der letzte Strahl
#define LIGHTSRC_ERROR -1        ///< Fehler ist aufgetreten
#define Z0 376.730313461         ///< Wellenwiderstand des freien Raumes 

#define LIGHTSRC_POL_X 0
#define LIGHTSRC_POL_Y 1
#define LIGHTSRC_POL_Z 2
#define LIGHTSRC_POL_USER_DEFINED 0 



/*!
LightSrc: Allgemeine Klasse zur Darstellung einer Lichtquelle
Die Lichtquelle befindet sich an der Position Pos. Sie wird durch eine quadratische Fläche repräsentiert,
von der die Strahlen loslaufen. N gibt dabei die Anzahl von Strahlen pro Seitenlänge der Fläche an, d.h.
es werden N^2 Strahlen berechnet. Die Definition der Strahlrichtung wird in den abgleiteten Klassen spezifiziert.
*/
class LightSrc
	
{
public:	 
  void reset();  ///< Alles zurücksetzen (Zählung beginnt von vorne)
  ~LightSrc(void);
  LightSrc() { type = LIGHTSRC_ERROR; Ein = 0; AnzObjs = 0; P0 = 1.0; }
  LightSrc(const LightSrc &);
  void clearObjects(); ///< Objektliste löschen
  void ObjectList (int Anz, Form **Ein);  ///< Objektliste setzen
  virtual int next(Ray_pow &S)=0; ///< Nächsten Strahl initialisieren (IStrahl-Klasse)
  virtual int next(IStrahl &S)=0; ///< Nächsten Strahl initialisieren (IStrahl-Klasse)
  virtual int next(Ray &S)=0;  ///< Nächsten Strahl initialisieren (Strahl-Klasse)
  void binRead (ifstream &is); 
  void binWrite (ofstream &os); ///< binäres Schreiben in die Datei os
  virtual void binWriteItem (ofstream &os) = 0; ///< binäres Schreiben in die Datei os (wird in den entsprechenden abgeleiteten Typen spezifiziert)
  virtual void binReadItem (ifstream &os) = 0; ///< binäres Lesen aus der Datei is (wird in den entsprechenden abgeleiteten Typen spezifiziert)
  int getNumObjs() { return AnzObjs; } ///< Anzahl Objekte ausgeben
  Form *getObject(int i) { if ((i<0) || (i>AnzObjs)) return NULL; return Ein[i]; }  ///< gibt Objekt zurück
  void setObject(Form *O, int i=-1); ///< ändert Objekt i in der Objektliste 
  int rayType(){return raytype; }  ///< gibt den Strahltyp zurück
  void setR0(double r0){ this->r0 = r0; }
  double getDensity () { return density; }
  void setD(double D) {
	this->density=D/((double) N);
	this->D=D;
	reset();}
  int getAnzRays() { return N; } ///< gibt Anzahl N der Strahlen zurück
  void setAnzRays(int N) ///< setzt die Anzahl N der Strahlen 
{
   this->N=N;
   density = D / ((double)N);
   cout << " N=" << this->N << endl;
   reset();
}
  void setPos(Vector<double> P); ///< setzt die Position der Lichtquelle
  void setN0(complex<double> n0){ this->n0=n0;} //
  Form **Ein;         ///< Objekte
  Vector<double> getk() {return k;}
// protected :
  	Vector<double> Pos; ///< Position der Lichtquelle
	int type;           ///< Typ der Lichtquelle
	 double P0;        ///< Leistung
    double density;     ///< Strahlen pro Raumrichtung = Strahldichte (=Abstand zwischen 2 Strahlen)
    Vector<double> k;   ///< Strahlrichtung;    
	int N,i1,i2;          ///< Indizes des zuletzt gestarteten Strahls (-1 -> noch kein Strahl gestartet)	
	Vector<complex<double> > Pol; ///< Polarisation der Quelle
	double r0;          ///< Radius der umgebenden "Weltkugel" (=Berechnungsraum)
	double wvl;         ///< Wellenlänge
	int AnzObjs;        ///< Anzahl Objekte
	complex<double> n0; ///< Brechungsindex Umgebung
  
	double D;           ///< Breite des Strahlfeldes am Anfang
	Vector<double> e1,e2;  ///< Einheitsvektoren, die den Bereich aufspannen, von denen die Strahlen loslaufen
	int raytype;        ///< Strahltyp : ray oder ISTRAHL (=RAY oder IRAY)
	int polType;        ///< Polarisationsrichtung (s.o.)  
friend class LightSrcPlane; 
friend class LightSrcGauss;
friend ostream& operator << (ostream &os, LightSrc *ls);
};



/*!
Klasse zur Darstellung einer ebene Welle als Lichtquelle
*/
class LightSrcPlane : public LightSrc
	
{
public:
	LightSrcPlane(void );	
	///< Konstruktor: mit folgenden Paramtern:
	///< Vector<double> Pos : Position der Lichtquelle 
	///< Vector<double> k   : Richtung der Lichtquelle
	///< int N : Anzahl der Strahlen in 
	LightSrcPlane (Vector<double> Pos, int N, double wvl, double r0=100.0, double D=100.0,  Vector<complex<double> > Pol=Vector<complex<double> >(0.0,1.0,0.0), int raytype=LIGHTSRC_RAYTYPE_IRAY);
	 LightSrcPlane(const LightSrcPlane &);
	~LightSrcPlane(void){};
	
	int next(IStrahl &S);
	int next(Ray &S);
	int next(Ray_pow &S);
	void binWriteItem (ofstream &os) { /* to be implemented !!! */ }
    void binReadItem (ifstream &os) { /* to be implemented !!! */ }
	
    // void turnSrc // to be done !!!

protected :
};

/**
 Klasse zur Beschreibung einer Lichquelle. Die Richtung des fokussierten Strahls ergibt sich aus der Quellposition und der Fokusposition
 Die elektrische Feldstärke wird berechnet gemäss:
 \f$\vec E(r,z)=\vec E_0 \frac{w_0}{w(z)}\cdot e^{r^2/w^2(z)}\cdot e^{-ik\frac{r^2}{2R(z)}}\cdot e^{i(\zeta(z)-kz)}\f$
*/
class LightSrcGauss : public LightSrc
{
public:
	void binWriteItem (ofstream &os);
    void binReadItem (ifstream &is);
	int next(Ray_pow &S);
	int next(IStrahl &S);	
	int next(Ray &S);	
	void initGauss(Vector<double> &P, Vector<double> &k, Vector<complex<double> > &E); ///< Initialisierung eines Teilstrahls
	LightSrcGauss(void);
//	LightSrcGauss(const LightSrc &);
	LightSrcGauss(const LightSrcGauss &);
    LightSrcGauss (Vector<double> Pos,  int N, double wvl, double w0, Vector<double> focuspos, double D=1.0,  Vector<complex<double> > Pol=Vector<complex<double> >(0.0,1.0,0.0), int raytype=LIGHTSRC_RAYTYPE_IRAY, double r0=1.0);
	void setW0(double w0) {this->w0=w0; calcz0(); reset();}
	void setFocuspos(Vector<double> fp) {focuspos=fp; f=abs(Pos-focuspos); reset(); }
	void setNA(double na); ///< setzt NA und berechnet D und w0 sowie z0 neu
	void setWvl(double wvl); ///< setzt (Vakuum-)Wellenlänge
	Vector<double> getFocuspos() {return focuspos; } ///< gibt Fokusposition 
	void setk(Vector<double> k) {this->k=k; reset(); }	
	double calcz0() { z0=M_PI*w0*w0/wvl; return z0; }
	void reset() 
	{
	 i1=0;
	 i2=0;
	 calcNormfak();
	}
	double calcw(double z)
	// Berechnung der Strahltaille, z: Abstand zum Fokuspunkt
	// Gibt den Wert zurück und speichert in gleichzeitig in die lokale Variable w
	{
      w=w0*sqrt(1.0+z*z/(z0*z0));
	  return w;
	}

	void calcNormfak()
    // Berechnung des Normierung
	{
	 double l=abs(Pos-focuspos);
	 calcz0();	 
	 calcw(l);
	 Normfak=Z0*sqrt(2.0/M_PI)/(n0*l*M_PI)*w/w0/w0*exp(-2.0*l*l/w/w);
	 Normfak=Normfak/((double) N*N)*l*l;
	}
//protected : 
	double w0;        ///< Taillendurchmesser (fiktiv !)
   
	double f;         ///< Abstand Startebene <-> Fokus
	double getNA() { return NA; }
protected :
	complex<double> Normfak;
	Vector<double> focuspos; ///< Ort des Fokalpunktes	
	double z0;     
	Vector<double> k; ///< Richtung des Gaußstrahls
	double w;
	double NA;    ///< Numerische Apertur (normiert auf Brechungsindex)
};

void binWriteLSList (ofstream &os, int nLS, LightSrc **ls);
void binReadLSList (ifstream &is, int nLS, LightSrc **&ls);
void copyLightSrcList (LightSrc **&d, LightSrc **s, int anz);

#endif