#include "surface.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>

//#define TEST_CULL

dreieck::dreieck()
{ 
//  P1=Vector<double>(0.0,0.0,0.0);
//  P2=Vector<double>(0.0,0.0,0.0);
//  P3=Vector<double>(0.0,0.0,0.0);
//  u=0.0;
//  v=0.0;
}

dreieck::dreieck(Vector<double> ip1, Vector<double> ip2, Vector<double> ip3)
{
  P[0] = ip1; P[1] = ip2; P[2] = ip3;
  setnorm();
}

void dreieck::setnorm()
{
  double hilf;
  n=(P[0]-P[2])%(P[1]-P[2]);
  hilf = P[2]*n;
#ifdef __GNUC__
  n=n*copysign(1.0,hilf);
#else
  n=n*_copysign(1.0,hilf);
#endif
  n/=abs(n);
}

Vector<double> dreieck::getnorm(void)
{
// setnorm();
 return n;
}

dreieck::dreieck(Vector<double> ip1, Vector<double> ip2, Vector<double> ip3, Vector<double> P0)
{  
  P[0] = ip1+P0; P[1] = ip2+P0; P[2] = ip3+P0;
  u=0.0; v=0.0;
  setnorm();
}

int dreieck::schnittpunkt(Vector<double> r, Vector<double> k, Vector<double> &p, double eps)
{
  Vector<double> edge1, edge2, pvec, qvec, tvec;
  double det, t, inv_det;

  edge1 = P[1] - P[0]; edge2 = P[2] - P[0];

  pvec = k%edge2;
  det = edge1*pvec;

#ifdef TEST_CULL

  if (det < eps)
    return 0;

  tvec = r - P[0];
  u = tvec*pvec;

  if (u < 0.0 || u > det)
    return 0;

  qvec = tvec%edge1;
  v = k*qvec;

  if (v < 0.0 || u+v > det)
    return 0;

  t=edge2*qvec;

  inv_det=1.0/det;

  t*=inv_det;
  u*=inv_det;
  v*=inv_det;

#else

  if (det > -eps  && det < eps)
    return 0;

  inv_det = 1.0/det;
  tvec = r - P[0];

  u = tvec*pvec*inv_det;

  if (u < 0.0 || u > 1.0)
    return 0;

  qvec = tvec%edge1;

  v = k*qvec*inv_det;

  if (v<0.0 || u+v > 1.0)
    return 0;

  t = edge2*qvec*inv_det;

#endif


  p = (1-u-v)*P[0] + u*P[1] + v*P[2];

  if (t < eps)
   return 0;

  return 1;
}

int dreieck::schnittpunkt(Vector<double> r, Vector<double> k,
                          double &t, Vector<double> &p,double eps)
{
  Vector<double> edge1, edge2, pvec, qvec, tvec;
  double det, inv_det;
 
  edge1 = P[1] - P[0]; edge2 = P[2] - P[0];
 
  pvec = k%edge2;
  det = edge1*pvec;
 
#ifdef TEST_CULL

  if (det < eps)
    return 0;
 
  tvec = r - P[0];
  u = tvec*pvec;
 
  if (u < 0.0 || u > det)
    return 0;

  qvec = tvec%edge1;
  v = k*qvec;

  if (v < 0.0 || u+v > det)
    return 0;

  t=edge2*qvec;

  inv_det=1.0/det;

  t*=inv_det;
  u*=inv_det;
  v*=inv_det;

#else

 if (det > -eps  && det < eps)
    return 0;

  inv_det = 1.0/det;
  tvec = r - P[0];

  u = tvec*pvec*inv_det;

  if (u < 0.0 || u > 1.0)
    return 0;

  qvec = tvec%edge1;

  v = k*qvec*inv_det;

  if (v<0.0 || u+v > 1.0)
    return 0;

  t = edge2*qvec*inv_det;

#endif


  p = (1-u-v)*P[0] + u*P[1] + v*P[2];


  if (t < eps)
  return 0;

  return 1;
}

dreieck::~dreieck()
{
}

Vector<double>& dreieck::operator[](int i)
{
 return P[i];
}

const Vector<double>& dreieck::operator[](int i)
const
{
 return P[i];
}

void dreieck::binWrite (ofstream &os)
{
 for (int i=0; i<3; i++)
  P[i].binWrite(os);
}

void dreieck::binRead (ifstream &is)
{
 for (int i=0; i<3; i++)
  P[i].binRead(is);
}



ostream& operator << (ostream &os, const dreieck &dr)
{
 os << dr[0] << "  " << dr[1] << "  " << dr[2]; // << "  " << dr.n;
 return os;
}

surface::surface(int anz)
{
  FName=0;
  r0=1.0;
  sf=1.0;
  anzp=anz;
  S=new dreieck[anzp];
  P=zero;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
  type=SURFACE;
}

surface::surface()
{
  FName=0;
  r0=1.0;
  sf=1.0;
  anzp=0;
  S=NULL;
  P = zero;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
  type=SURFACE;
}

surface::surface (const Form &F)
{
 r0=1.0;
 sf=1.0;
 type=SURFACE;
 FName=0;
}

surface::surface(const surface &Su):Form(Su)
{
  
  anzp=Su.anzp;
  sf=Su.sf;
  r0=Su.r0;
  sf=1.0;  
  P = zero;
  currentnorm=Su.currentnorm;
  currentIndex = Su.currentIndex;
  S=new dreieck[anzp];
  P=Su.P;
  for (int i=0;i<anzp;i++)
   S[i]=Su.S[i];
  type=SURFACE;  
//  if (FName!=0) delete[] FName;  
  FName=new char[strlen(Su.FName)+1];
  memcpy (FName,Su.FName,strlen(Su.FName)+1);
//   FName="HALLO";  
   
//     sprintf (FName,"%s",Su.FName);
}

surface::surface(Vector<double> Oh)
{
  FName=0; 
  anzp=0;
  S=NULL;
  P = Oh;
  r0=1.0;
  type=SURFACE;
}


surface::surface(const Vector<double> &P,
                complex<double>  n,
          const Matrix<complex<double> > alpha,
          const Vector<double> &Ex,
          const Vector<double> &Ey,
          const Vector<double> &Ez
         )
: Form (P, n, alpha, Ex, Ey, Ez, SURFACE)
{
  FName=0;
 r0=1.0;
  sf=1.0;
  anzp=0;
  S=NULL;
  currentnorm=Vector<double>(0,0,0);
  currentIndex = -1;
 trafo (Ex,Ey,Ez,H,R);
 this->P=P;
 type=SURFACE;
} 

surface::surface(const Vector<double> &P,
                complex<double>  n,
          int anz, dreieck* list,
          const Matrix<complex<double> > alpha,
          const Vector<double> &Ex,
          const Vector<double> &Ey,
          const Vector<double> &Ez)
: Form (P, n, alpha, Ex, Ey, Ez,SURFACE)
{
 FName=0;
 trafo (Ex,Ey,Ez,H,R);
 this->P=P;
 this->anzp=anz;
 S=new dreieck[anzp];
 for(int i=0;i<anzp;i++)
 this->S[i]=list[i];
 type=SURFACE;
 initQuad();
}

surface::~surface()
{
 clearS();
}


int surface::createsurface()
{
  double P1x, P1y, P1z, P2x, P2y, P2z, P3x, P3y, P3z;
  Vector<double> P1, P2, P3;
  cout  << "Anzahl Flaechen:";
  cin >> anzp;
  delete[] S;
  S = new dreieck[anzp];

  for(int i=0;i<anzp;i++){
    cout << "Flaeche" << i+1 << endl;
    cout << "P1x:" << endl;
    cin >> P1x;
    cout << "P1y:" << endl;
    cin >> P1y;
    cout << "P1z:" << endl;
    cin >> P1z;
    cout << "P2x:" << endl;
    cin >> P2x;
    cout << "P2y:" << endl;
    cin >> P2y;
    cout << "P2z:" << endl;
    cin >> P2z;
    cout << "P3x:" << endl;
    cin >> P3x;
    cout << "P3y:" << endl;
    cin >> P3y;
    cout << "P3z:" << endl;
    cin >> P3z;
    P1 = Vector<double>(P1x,P1y,P1z);
    P2 = Vector<double>(P2x,P2y,P2z);
    P3 = Vector<double>(P3x,P3y,P3z);
    S[i]=dreieck(P1,P2,P3);
  }
  return 0;
}

int surface::createsurface(char* FName)
{
  Vector<double> P1, P2, P3;

  ifstream is;
  if (this->FName!=0) delete this->FName;
  this->FName=new char[strlen(FName)+1];
  strcpy (this->FName,FName);
  is.open(FName);
  is >> anzp;
 // cout << "Create anzp:" << anzp << endl;

  S = new dreieck[anzp];

  for(int i=0;i<anzp;i++)
  {
    is >> P1 >> P2 >> P3;
    S[i] = dreieck(P1*r0,P2*r0,P3*r0); 
  //  cout << "r0:" << r0 << ", S[" << i <<"]:" << S[i] << endl;
  }
//  cout << "anzp nach:" << anzp << endl;

  initQuad();
   setCenter2CoM();
  return 0;
}

int surface::importBinSTL(char *FName)
{
	ifstream is;
	int dummy;
	int i,j;
	int anz;
	float data;
	char str[255];
	Vector<double> P1,P2,P3,n;
	Vector<double> cm;
	cout << "---------------------------- IMPORT STL-FILE --------------------------------" << endl;
	if (anzp!=0) delete[] S;

	is.open (FName,ios::binary);
	if (is.good())
	{
	is.read (str,80);
	anz=readLE_int32(is);
	if (this->FName!=0) delete[] this->FName;
	this->FName=new char[strlen(FName)+1];
	// sprintf (this->FName,"%s",FName);
	strcpy (this->FName,FName);
	S=new dreieck[anz];
	anzp=anz;
	cout << "Lese " << anz << "  Dreiecke" << endl;

	for (i=0; i<anz && !is.eof(); i++)
	{
		// is.read(str,12);  // Dreiecksnormale --> wird eh neu berechnet 

		for (j=0; j<3; j++)   // Punkt 1
		{  
         n[j]=readLE_float32(is);
		}
		
		for (j=0; j<3; j++)   // Punkt 1
		{  
         P1[j]=readLE_float32(is);
		 cm[j]+=P1[j];
		}

		for (j=0; j<3; j++)   // Punkt 2
		{
		 P2[j]=readLE_float32(is);
		 cm[j]+=P2[j];
		}

		for (j=0; j<3; j++)   // Punkt 3
		{
		 P3[j]=readLE_float32(is);
		 cm[j]+=P3[j];
		}

		cm=cm/anz/3.0;
		is.read(str,2);
//		cout << P1 << "   " << P2 << "   " << P3 <<  "    " << n << endl;

//		S[i]=dreieck(P1-cm,P2-cm,P3-cm);
		S[i]=dreieck(P1,P2,P3);
		S[i].setnorm(n); // Glauben wir mal, dass die Oberflächennormale im STL-File richtig ist !
	}
	initQuad();
	is.close();
	setCenter2CoM();
	cout <<  "------------------------------- IMPORT ENDE ---------------------------------" << endl;
	return 0;
	}
	else return -1;
}

bool surface::next(const Vector<double> &r, const Vector<double> &k, Vector<double> &p, const int inside)
{
// Angenommene Koordinatensysteme:
//
// r     : Globales Koordinatensystem
// S[i].P: Einschlusskoordinatensystem
// rhilf : r im Einschlusskoordinatensystem
// khilf : k im Einschlusskoordinatensystem
// p     : Globales Koordinatensystem
//
  double thilf,t;
  Vector<double> rhilf, philf, khilf;
  bool found=false;
  rhilf = H*(r - P);
  khilf = H*k;
  currentIndex = -1;
  for(int i=0;i<anzp;i++)
  {
    if(S[i].schnittpunkt(rhilf,khilf,t,p))
    {
      // found=true;
       if((!found) && (t>0))
       {
        thilf=t; philf=p; currentnorm=S[i].n;
	currentIndex = i;
        found=true;
       }
       else
       {
        if((t<thilf)&&(t>1.e-8))
        {
           thilf=t;
           philf=p;
           currentnorm=S[i].n;
	   currentIndex = i;
        }
       }
    }
//    else
//    {
//      cout << i << endl;
//      cout << "nix\n";
//    }
  }
  if (found)
  {
   p=R*philf+P;
 //  cout << "next->currentIndex=" << currentIndex << endl;
   return found;
  }
  else
  { 
//   cout <<"Kein Schnittpunkt!" << endl;
//   p=Vector<double> (-1E20,-1E20,-1E20);
   return found;
  }

}

bool surface::isInside(const Vector<double> &p)
{
  double hilf;
  Vector<double>  hilf2;
  cout << "XXX p:" << p << endl;
  for(int i=0;i<anzp;i++)
  {
    hilf2 = S[i].P[2]-p;
    hilf = (S[i].n*hilf2);
    cout << "Dreieck" << i <<": " << S[i] << endl;
    cout << "n:" << S[i].n << endl;
    cout << "hilf2:" << hilf2 << endl;
    cout << "hilf:" << hilf << endl;
    if(hilf<0.0)
      return false;
  }
  return true;
}

void surface::initQuad()
{
  double xmin,ymin,zmin,xmax,ymax,zmax;
  dreieck Sh;

  Sh = R*S[0]+P;
  xmin=Sh.P[0][0]; ymin = Sh.P[0][1]; zmin=Sh.P[0][2];
  xmax=Sh.P[0][0]; ymax = Sh.P[0][1]; zmax=Sh.P[0][2];

  for(int i=0;i<anzp;i++)
  {
   Sh = R*S[i]+P;
   if(Sh.P[0][0]<xmin)
     xmin = Sh.P[0][0];
   else if(Sh.P[0][0]>xmax)
     xmax = Sh.P[0][0];
   if(Sh.P[0][1]<ymin)
     ymin = Sh.P[0][1];
   else if(Sh.P[0][1]>ymax)
     ymax = Sh.P[0][1];
   if(Sh.P[0][2]<zmin)
     zmin = Sh.P[0][2];
   else if(Sh.P[0][2]>zmax)
     zmax = Sh.P[0][2];

   if(Sh.P[1][0]<xmin)
     xmin = Sh.P[1][0];
   else if(Sh.P[1][0]>xmax)
     xmax = Sh.P[1][0];
   if(Sh.P[1][1]<ymin)
     ymin = Sh.P[1][1];
   else if(Sh.P[1][1]>ymax)
     ymax = Sh.P[1][1];
   if(Sh.P[1][2]<zmin)
     zmin = Sh.P[1][2];
   else if(Sh.P[1][2]>zmax)
     zmax = Sh.P[1][2];

   if(Sh.P[2][0]<xmin)
     xmin = Sh.P[2][0];
   else if(Sh.P[2][0]>xmax)
     xmax = Sh.P[2][0];
   if(Sh.P[2][1]<ymin)
     ymin = Sh.P[2][1];
   else if(Sh.P[2][1]>ymax)
     ymax = Sh.P[2][1];
   if(Sh.P[2][2]<zmin)
     zmin = Sh.P[2][2];
   else if(Sh.P[2][2]>zmax)
     zmax = Sh.P[2][2];
  }

 /* cout << "xmin=" << xmin << "   xmax=" << xmax;
  cout << "  ymin=" << ymin << "   ymax=" << ymax;
  cout << "  zmin=" << zmin << "   zmax=" << zmax << endl; */
  pul=Vector<double>(xmin,ymin,zmin)+P;
  por=Vector<double>(xmax,ymax,zmax)+P;
 // cout << "pul=" << pul << "    por=" << por << endl;
}

Vector<double> surface::norm (const Vector<double> &dummy)
{
 // Dummy Argument ist notwendig aufgrund des Designs der Formklasse
 return R*currentnorm;
}

void surface::setr0(double rneu)
{
  cout << "surface::setr0" << endl;
  for(int i=0;i<anzp;i++)
  { 
    S[i].P[0]=S[i].P[0]*rneu/r0;
    S[i].P[1]=S[i].P[1]*rneu/r0;
    S[i].P[2]=S[i].P[2]*rneu/r0;
  }
  r0=rneu;
  initQuad();
}

surface surface::nosurface()
{
  surface h;
  h=*this;
  h.clearS();
  return h; 
}

void surface::clearS()
{
  if (anzp!=0)
  {
    delete[] S;
    S = NULL;
    anzp=0;
  }
}

void surface::addS(dreieck* list,int anz)
{
  cout << "addieren" << endl;
  if(anzp!=0)
   clearS();

  anzp=anz;
  
  S = new dreieck[anzp];
  
  for(int i=0;i<anzp;i++)
  {
    S[i]=list[i];
    cout << S[i] << endl;
  }

}



ostream& operator << (ostream &os, const surface &su)
{
	os << "Pos=" << su.P << endl;
	os << "Winkel:" << su.Ealpha << "," << su.Ebeta << "," << su.Egamma << endl; 
 os << "anzp:" << su.anzp << endl;
 if(su.anzp>0)
 {
  os << "[";
  for (int i=0;i<su.anzp-1;i++)
    os  << su.S[i] << "\n";
  os << su.S[su.anzp-1] << "]";
 }
 else
 {
  os <<"[]";
 } 
 return os;
}

dreieck operator + (const dreieck &dr, const Vector<double> &v)
{
 dreieck h;
 h[0] = dr[0]+v;
 h[1] = dr[1]+v; 
 h[2] = dr[2]+v;
 h.setnorm();

 return h;
}

dreieck operator + (const Vector<double> &v, const dreieck &dr)
{
 dreieck h;
 h[0] = v+dr[0];
 h[1] = v+dr[1];
 h[2] = v+dr[2];
 h.setnorm();


 return h;
}


dreieck operator - (const dreieck &dr, const Vector<double> &v)
{
 dreieck h;
 h[0] = dr[0]-v;
 h[1] = dr[1]-v;
 h[2] = dr[2]-v;
 h.setnorm();

 return h;
}

dreieck operator - (const Vector<double> &v, const dreieck &dr)
{
 dreieck h;
 h[0] = v-dr[0];
 h[1] = v-dr[1];
 h[2] = v-dr[2];
 h.setnorm();

 return h;
}

dreieck operator * (const Matrix<double> &M, const dreieck &dr)
{
 dreieck h;
 h[0] = M*dr[0];
 h[1] = M*dr[1];
 h[2] = M*dr[2];
 h.setnorm();
 
 return h;
}

dreieck operator / (const dreieck &dr, double a)
{
 dreieck h;
 h[0] = dr[0]/a;
 h[1] = dr[1]/a;
 h[2] = dr[2]/a;
 h.setnorm();
 return h;
}


dreieck operator * (const dreieck &dr, double a)
{
 dreieck h;
 h[0] = a*dr[0];
 h[1] = a*dr[1];
 h[2] = a*dr[2];
 h.setnorm();
 return h;
}

dreieck operator * (double a, const dreieck &dr)
{
 dreieck h;
 h[0] = a*dr[0];
 h[1] = a*dr[1];
 h[2] = a*dr[2];
 h.setnorm();
 return h;
}

dreieck& dreieck::operator = (const dreieck &dr)
{
 this->P[0][0]=dr.P[0][0];
 this->P[0][1]=dr.P[0][1];
 this->P[0][2]=dr.P[0][2];
 this->P[1][0]=dr.P[1][0];
 this->P[1][1]=dr.P[1][1];
 this->P[1][2]=dr.P[1][2];
 this->P[2][0]=dr.P[2][0];
 this->P[2][1]=dr.P[2][1];
 this->P[2][2]=dr.P[2][2];
 this->n[0]=dr.n[0];
 this->n[1]=dr.n[1];
 this->n[2]=dr.n[2];
 return *this;
}

surface& surface::operator = (const surface &s)
{
  this->P=s.P;
  this->H=s.H;
  this->R=s.R;
  this->n=s.n;
  this->alpha=s.alpha;
  this->type=s.type;
  this->pul=s.pul;
  this->por=s.por;
  for(int i=0;i<3;i++)
  this->e[i]=s.e[i];
  this->Ealpha=s.Ealpha;
  this->Ebeta=s.Ebeta;
  this->Egamma=s.Egamma;
  this->r0=s.r0;
  this->anzp=s.anzp;
  this->S=new dreieck[this->anzp];
  for(int i=0;i<this->anzp;i++)
  this->S[i]=s.S[i];
  this->sf=s.sf;
  this->currentnorm=s.currentnorm;
  this->currentIndex = s.currentIndex;
   sprintf (this->FName,"%s",FName);
  //strcpy (this->FName,FName);
  return *this;
}

surface operator + (const surface &s, const Vector<double> &v)
{

 surface h(s);

 for (int i=0;i<s.anzp;i++)
 {
   h.S[i]=s.S[i]+v;
 } 
 return h;
}

surface operator - (const surface &s, const Vector<double> &v)
{

 surface h(s);
 
 for (int i=0;i<s.anzp;i++)
 {
   h.S[i]=s.S[i]-v;
 } 
 return h;
}

surface operator * (const Matrix<double> &M, const surface &s)
{
 surface h(s);

 for (int i=0;i<s.anzp;i++)
 {
   h.S[i]=M*s.S[i];
 }
 return h;
}
/** No descriptions */
void surface::scale (double sf)
{
 for (int i=0; i<anzp; i++)
  S[i]=S[i]*sf/this->sf;
 this->sf=sf;
}
/** No descriptions */
double surface::isInHost(void)
{
 double D,Pk,l;
 double erg=1;
 double r=0;
 Vector<double> hv,k,ergv;
 if (abs(P)>1.0) return -1;
 for (int i=0; i<anzp; i++)
 {
  for (int j=0; j<3; j++)
  {
   hv=R*S[i][j]+P;
   r=abs(hv);
   if ((r>erg) && (r>1.0)) { erg=r; ergv=hv; }
  }
 }
 if (erg>1.0)
 {
  k=hv-P;
  Pk=P*k;
  D=Pk*Pk-abs2(k)*(abs2(P)-1.0);
  l=-Pk+sqrt(D)/(abs2(k));
  return l/abs(ergv);
 }
 return 1.0/erg;
}

void surface::binWrite (ofstream &os)
{
 size_t strl;
 P.binWrite(os);
 H.binWrite(os);
 R.binWrite(os);
 os.write ((char *) &n, (char) sizeof (n)); 
 alpha.binWrite(os);
 pul.binWrite (os);
 por.binWrite (os);
 for (int i=0; i<3; i++)
  e[i].binWrite(os);
 os.write ((char *) &Ealpha,(char) sizeof (Ealpha));
 os.write ((char *) &Ebeta,(char) sizeof (Ebeta));
 os.write ((char *) &Egamma,(char) sizeof (Egamma));
 os.write ((char *) &sf,(char) sizeof (sf));
 os.write ((char *) &r0,(char) sizeof (r0));
 os.write ((char *) &anzp, (char) sizeof(anzp));
 strl=strlen(FName);
 os.write ((char *) &strl, sizeof(strl));
 char c;
 for (int i=0; i<=strlen(FName); i++)
 {
  c=FName[i];
  cout << "c=" << c << endl;
 os.write ((char *) &c, 1);
 }
 for (int i=0; i<anzp; i++)
  S[i].binWrite(os);  
}

void surface::binRead (ifstream &is)
{
 type=SURFACE;
 P.binRead(is);
 H.binRead(is);
 R.binRead(is);
 is.read((char *) &n, (char) sizeof(n)); 
 alpha.binRead(is);
 pul.binRead (is);
 por.binRead (is);
 for (int i=0; i<3; i++)
  e[i].binRead(is);
 is.read ((char *) &Ealpha,(char) sizeof (Ealpha));
 is.read ((char *) &Ebeta,(char) sizeof (Ebeta));
 is.read ((char *) &Egamma,(char) sizeof (Egamma));
 is.read ((char *) &sf,(char) sizeof (sf));
 is.read ((char *) &r0,(char) sizeof (r0));
 clearS();
 is.read ((char *) &anzp,(char) sizeof (anzp));
 size_t strl;
 is.read ((char *) &strl, (char) sizeof(strl));
 if (FName!=0) delete FName;
 FName=new char[strl+1];
 char c;
 for (int i=0; i<=strl; i++)
 {
  is.read ((char *) &c, 1);
  FName[i]=c;
 } 
//  cout << "FName=" << FName << endl; 
 S=new dreieck[anzp];
 for (int i=0; i<anzp; i++)
 {
  S[i].binRead(is);
  S[i].setnorm();
 } 
}


void surface::exportSRF (char *FName)
{
 ofstream os;
 os.open (FName);
 os << anzp << endl;
 for (int i=0; i<anzp; i++)
  os << S[i]/r0 << endl;
 os.close();
}

/*!
    \fn surface::calcVolume()
 */
double surface::Volume() {return calcVolume(); }

double surface::calcVolume()
{
 /* double F,Vg,absn,h;
  Vector<double> n;
  
  if (anzp==0) return 0.0; 
  Vg=0.0;
  for (int i=0; i<anzp; i++)
  {
   n=(S[i][1]-S[i][0]) % (S[i][2]-S[i][0]);
   absn=abs(n);
   n/=absn;
   F=absn/2.0;
   h=fabs((S[i][0]-P)*n);   
   Vg+=F*h/3.0;
//    cout << "h=" << h << endl;
  }
  return Vg;*/
  double Vg=0;
  double A;
  Vector<double> vsum;
  
  for (int i=0; i<anzp; i++)
  {
    A=S[i].area();
	vsum=S[i].P[0]+S[i].P[1]+S[i].P[2];
	Vg+=A/3.0*S[i].getnorm()*vsum;	
  }
  return Vg/3.0;
}
void surface::setP (Vector<double> r)
{
 //Vector<double> dP=P-r;
 P=r;
 /*for (int i=0; i<anzp; i++)
   for (int j=0; j<3; j++)
  S[i].P[j]+=dP;*/ 
 initQuad();
}


Vector<double> surface::calcCoM()
{
	Vector<double> n,Fx,Fy,Fz;
	double int_x=0,int_y=0,int_z=0;
	double F,x,y,z,V=calcVolume(),Fges=0;
	for (int i=0; i<anzp; i++)
	{
		F=S[i].area();
		Fges+=F;
		x=S[i][0][0];
		y=S[i][0][1];
		z=S[i][0][2];
		n=S[i].getnorm();
		Fx=Vector<double> (0.5*x*x,0,0);
		Fy=Vector<double> (0,0.5*y*y,0);
		Fz=Vector<double> (0,0,0.5*z*z);
		int_x+=F*Fx*n;
		int_y+=F*Fy*n;
		int_z+=F*Fz*n;
	}
	cout << "Gesamtfläche :" << Fges << endl;
	Vector<double> P=Vector<double>(int_x,int_y,int_z)/V;
	return P;
}

void surface::setCenter(Vector<double> P)
{
	cout << "P= " << P << "     anzp=" << anzp << endl;
	for (int i=0; i<anzp; i++)
	{
		S[i].P[0]=S[i].P[0]-P;
		S[i].P[1]=S[i].P[1]-P;
		S[i].P[2]=S[i].P[2]-P;
	}
}

inline void Subexpressions (double w0,double w1, double w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2)
{
	double temp0=w0+w1; f1=temp0+w2;
	double temp1=w0*w0;
	double temp2=temp1 + w1*temp0;
	f2=temp2+w2*f1; f3=w0*temp1+w1*temp2+w2*f2;
	g0=f2+w0*(f1+w0); g1=f2+w1*(f1+w1); g2=f2+w2*(f1+w2);
}

Matrix<double> surface::computeInertia()
{
	double f1x,f2x,f3x,g0x,g1x,g2x;
	double f1y,f2y,f3y,g0y,g1y,g2y;
	double f1z,f2z,f3z,g0z,g1z,g2z;
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	Vector<double> d;
	Matrix <double> inertia;
		
	const double mult[10]={1.0/6.0,1.0/24.0,1.0/24.0,1.0/24.0,1.0/60.0,1.0/60.0,1.0/60.0,1.0/120.0,1.0/120.0,1.0/120.0};
	double intg[10]={0,0,0,0,0,0,0,0,0,0};
	for (int t=0; t<anzp; t++)
	{
		x0=S[t].P[0][0]; y0=S[t].P[0][1]; z0=S[t].P[0][2];
		x1=S[t].P[1][0]; y1=S[t].P[1][1]; z1=S[t].P[1][2];
		x2=S[t].P[2][0]; y2=S[t].P[2][1]; z2=S[t].P[2][2];

		d=(S[t].P[1]-S[t].P[0])%(S[t].P[2]-S[t].P[0]);
//		d=d*(S[t].n*d); // Sicherstellen, dass d auch nach außen zeigt !
		
		Subexpressions(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x);
		Subexpressions(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y);
		Subexpressions(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);

		intg[0]+=d[0]*f1x;
		intg[1]+=d[0]*f2x; intg[2]+=d[1]*f2y; intg[3]+=d[2]*f2z;
		intg[4]+=d[0]*f3x; intg[5]+=d[1]*f3y; intg[6]+=d[2]*f3z;
		intg[7]+=d[0]*(y0*g0x+y1*g1x+y2*g2x);
		intg[8]+=d[1]*(z0*g0y+z1*g1y+z2*g2y);
		intg[9]+=d[2]*(x0*g0z+x1*g1z+x2*g2z);
	}

	for (int i=0; i<10; i++) intg[i]*=mult[i];

	double mass=intg[0];

	// center of mass 
	Vector<double> cm=Vector<double> (intg[1],intg[2],intg[3])/mass;
	cout << "CoM=" << cm << endl;
	cout << "mass=" << mass << endl;
	inertia(0,0)=intg[5]+intg[6]-mass*(cm[1]*cm[1]+cm[2]*cm[2]);
	inertia(1,1)=intg[4]+intg[6]-mass*(cm[2]*cm[2]+cm[0]*cm[0]);
	inertia(2,2)=intg[4]+intg[5]-mass*(cm[0]*cm[0]+cm[1]*cm[1]);
	inertia(0,1)=-(intg[7]-mass*cm[0]*cm[1]);
	inertia(1,2)=-(intg[8]-mass*cm[0]*cm[2]);
	inertia(0,2)=-(intg[9]-mass*cm[2]*cm[0]);

	return inertia/mass;  
}

void surface::setCenter2CoM()
{
	Vector<double> CoM=calcCoM();
	setCenter(CoM);
}

surface generateHexagonCylinder(double a, double h, Matrix<double> M)
{
	int i;
	int N = 24; // 24 Dreiecke
	double s30 = 0.5; // sin 30°
	double c30 = SQRT3 / 2.0;
	dreieck *D = new dreieck[N];
	Vector<double> P[7],P1,P2,P3;
	Vector<double> hp = 0.5*h*ez;
	Vector<double> hm = -hp;

	P[0] = Vector<double>(0, a, 0);
	P[1] = Vector<double>(c30*a, s30*a, 0);
	P[2] = Vector<double>(c30*a, -s30*a, 0);
	P[3] = Vector<double>(0, -a, 0);
	P[4] = Vector<double>(-c30*a, -s30*a, 0);
	P[5] = Vector<double>(-c30*a, s30*a, 0);
	P[6] = P[0];
	for (i = 0; i < 6; i++)
	{
		D[i] = dreieck(M*(P[i] + hp), M*(P[i + 1] + hp), M*hp);  // oberer Deckel
                D[i].n = M*ez; 
		D[i + 6] = dreieck(M*(P[i] + hm),M*( P[i + 1] + hm), M*hm); // unterer Deckel
                D[i+6].n =-M*ez;
                
                P1=M*(P[i] + hp); 
                P2=M*(P[i] + hm); 
                P3=M*(P[i + 1] + hp);
		D[i + 12] = dreieck(P1, P2, P3 );
                D[i + 12].n=(P1-P2)%(P3-P2);
                D[i + 12].n/=abs(D[i + 12].n);

                P1=M*(P[i + 1] + hp);
                P2=M*(P[i + 1] + hm); 
                P3=M*(P[i] + hm);
		D[i + 18] = dreieck(P1,P2,P3);
                D[i + 18].n=-(P1-P2)%(P3-P2);
                D[i + 18].n/=abs(D[i + 18].n);
 
	}
	return 	surface(zero, 1.5, N, D);
}

