#include "funktion.h"
#include <sstream>
#include "funsurf.h"


Funsurf::Funsurf ()
: Form()
{
 type = FUNSURF;
 SurfGitter=0;
}

Funsurf::Funsurf (const Form &F) : Form (F)
{
 type = FUNSURF;
 SurfGitter=0;
}

Funsurf::Funsurf (const Funsurf &FS)
:Form (FS)
{
  type = FUNSURF;
  SurfGitter=0;
}

 Funsurf::Funsurf ( const string &S,
             const Vector<double> &P,
             const Matrix <complex<double> > &alpha,
             const double r0,
             const Vector <double> &Ex,
             const Vector <double> &Ey,
             const Vector <double> &Ez,
             const complex<double>  n)
:Form(P,n,alpha,Ex,Ey,Ez,FUNSURF)
{
 SurfGitter=0;
 Fkt=Funktion(S);
trafo(Ex,Ey,Ez,H,R);
this->P=P;
this->r0=r0;
type=FUNSURF;
}

Funsurf::Funsurf(string &S)
{
 SurfGitter=0;
 Fkt=Funktion(S);
type=FUNSURF;
}


Funsurf::Funsurf(
                          const Funktion &Fkt,
                          const Vector<double> &P,
                          const Matrix <complex<double> > &alpha,
                          const double r0,
                          const Vector <double> &Ex,
                          const Vector <double> &Ey,
                          const Vector <double> &Ez,
                          const complex<double>  n)
:Form(P,n,alpha,Ex,Ey,Ez,FUNSURF)
{
SurfGitter=0;
trafo(Ex,Ey,Ez,H,R);
this->P=P;
this->r0=r0;
this->Fkt=Fkt;
type=FUNSURF;
}




Funsurf::~Funsurf(){
}

bool Funsurf::next (const Vector<double> &Ps, const Vector<double> &K, Vector<double> &pout, const int inside)
{
int posl, posx, posy, posz;
string strx,stry,strz,strp[3],strk[3],x;
Vector<double> k,p=Ps-P;
Funktion f,df,fmaster;
k=K;
 p=H*p;
 k=H*K;
f=Fkt;
fmaster=Fkt;

stringstream ss;

ss << p[0];
ss >> strp[0];
ss.clear();

ss << k[0];
ss >> strk[0];
ss.clear();
// strx="("+strp[0]+"+l*"+"("+strk[0]+")"+")";
strx="("+strp[0]+"+"+strk[0]+"*l)";

ss << p[1];
ss >> strp[1];
ss.clear();
ss << k[1];
ss >> strk[1];
ss.clear();
//stry="("+strp[1]+"+l*"+"("+strk[1]+")"+")";
stry="("+strp[1]+"+"+strk[1]+"*l)";

ss << p[2];
ss >> strp[2];
ss.clear();
ss << k[2];
ss >> strk[2];
ss.clear();
//strz="("+strp[2]+"+l*"+"("+strk[2]+")"+")";
strz="("+strp[2]+"+"+strk[
2]+"*l)";

string str;
str=f.getfktm();

posx=0;
while(posx>=0)
  {
posx=str.find_first_of("x");
if (posx>=0) str=str.replace(posx,1,strx);
posx=str.find_first_of("x");
  }

posy=0;
while(posy>=0)
  {
posy=str.find_first_of("y");
 if (posy>=0) str=str.replace(posy,1,stry);
posy=str.find_first_of("y");
   }

posz=0;
while(posz>=0)
  {
posz=str.find_first_of("z");
 if (posz>=0)  str=str.replace(posz,1,strz);
posz=str.find_first_of("z");
  }

posl=0;
while (posl>=0)
 {
posl=str.find_first_of("l");
if(posl>=0) str=str.replace(posl,1,"x");
posl=str.find_first_of("l"); 
 }


f=Funktion(str);
/*
for(int i=0;i<=f.nS;i++)
{
cout<<"S"<<i<<"= "<<f.S[i]<<endl;
}*/
double ns;
bool nst;


// f.clean();
//nst=nullstelle(f,ns);
ns=f.Nullstelle(0,2*r0,r0/20.0,nst,1E-4*r0);

//df=f.partdiff("x");
pout[0]=p[0]+ns*k[0];
pout[1]=p[1]+ns*k[1];
pout[2]=p[2]+ns*k[2];


/*if (nst)
{cout<<"   Nullstelle bei: "<<pout<<"    abs: "<<abs(pout)<<"  ns= "<<ns<<" f(ns)= "<<f(ns)<<"    f=  "<<f<<endl;
cout<<"der strahl: p= "<<p<<"  abs(p)=  "<<abs(p)<<"   k=  "<<k<<endl;
cout<<"-----------------------------------------------------------"<<endl;   }
  */
return nst;
}

bool  Funsurf::nullstelle(Funktion f,double &ns)
{
double dx,y0,y1;
int n=3;
dx=r0/SW;
y0=f(dx);
y1=f(2*dx);
//cout << "r0=" << r0 << endl;
// cout << "f=" << f << endl;
while((y0*y1>0) && (n*dx<=2.0*r0))
{
 y0=y1;
 y1=f(n*dx);
// cout << "x=" << n*dx << "     y0=" << y0 << "     y1=" << y1 << endl;
 if((y0*y1)<=0)
 {
 ns=((n-0.5)*dx);
//  cout << "++++++:ns=" << ns << endl << "__________________" << endl;
 return true;
 }
 n=n+1;
}
// cout << "-----:ns=" << ns << endl << "__________________" << endl;
return false;
}

Vector<double> Funsurf::norm (const Vector<double> &p)
{
Vector<double> n;
Funktion f,h;
f=Fkt;

h=f.partdiff("x");
n[0]=h(p[0],p[1],p[2]);

h=f.partdiff("y");
n[1]=h(p[0],p[1],p[2]);

h=f.partdiff("z");
n[2]=h(p[0],p[1],p[2]);

n=n/abs(n);
return n;

}
/*
void Funsurf::scale(double sf)
{
stringstream ss;
string scf;
ss<<sf;
ss>>scf;
ss.clear();


} */
void Funsurf::BerechneSurfGitter(int ntheta, int nphi)
{
Vector<double> er,psurf;
string strp[3],strer[3],strx,stry,strz;
Funktion f;
double ns;
bool nst;
f=Fkt;
int posl, posx, posy, posz;
double dphi=2.0*M_PI/(double)(nphi-1);
double dtheta=M_PI/(double)(ntheta-1);
double theta,phi;
string fktm=f.getfktm();
stringstream ss;
if (SurfGitter!=0)
for (int i=0; i<this->ntheta; i++)
  delete[] SurfGitter[i];
delete SurfGitter;

this->ntheta=ntheta;
this->nphi=nphi;

SurfGitter=new Vector<double> *[ntheta];
for (int i=0; i<ntheta; i++)
  SurfGitter[i]=new Vector<double>[nphi];


for(int i=0 ; i<ntheta; i++)
 {
  theta=i*dtheta;
 er[2]=cos(theta);
 for(int j=0 ; j<nphi; j++)
  {
   phi=j*dphi;
  if(i==0 || i==ntheta)
   {
   er[0]= sin(theta);
   er[1]=0.0;
   }
  else
  {
  er[0]= cos(phi)*sin(theta);
  er[1]= sin(phi)*sin(theta);
  }



  ss.clear();
  ss << P[0];
  ss >> strp[0];
  ss.clear();

  ss << er[0];
  ss >> strer[0];
  ss.clear();
  strx="("+strp[0]+"+"+strer[0]+"*l)";

  ss << P[1];
  ss >> strp[1];
  ss.clear();
  ss << er[1];
  ss >> strer[1];
  ss.clear();
  stry="("+strp[1]+"+"+strer[1]+"*l)";

  ss << P[2];
  ss >> strp[2];
  ss.clear();
  ss << er[2];
  ss >> strer[2];
  ss.clear();
  strz="("+strp[2]+"+"+strer[2]+"*l)";

  string str;
  str=fktm;

posx=0;
while(posx>=0)
  {
posx=str.find_first_of("x");
if (posx>=0) str=str.replace(posx,1,strx);
posx=str.find_first_of("x");
  }

posy=0;
while(posy>=0)
  {
posy=str.find_first_of("y");
 if (posy>=0) str=str.replace(posy,1,stry);
posy=str.find_first_of("y");
   }

posz=0;
while(posz>=0)
  {
posz=str.find_first_of("z");
 if (posz>=0)  str=str.replace(posz,1,strz);
posz=str.find_first_of("z");
  }

posl=0;
while (posl>=0)
 {
posl=str.find_first_of("l");
if(posl>=0) str=str.replace(posl,1,"x");
posl=str.find_first_of("l");
 }


//cout << "str=" << str << endl;
f=Funktion(str);
str="";
//f.clean();
nst=nullstelle(f,ns);
psurf[0]=P[0]+ns*er[0];
psurf[1]=P[1]+ns*er[1];
psurf[2]=P[2]+ns*er[2];
//cout << "psurf=" << psurf <<endl;
SurfGitter[i][j]=psurf;

  }
 }   // end for
}


void Funsurf::binWrite (ofstream &os)
{
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
os.write ((char *) &nphi, (char) sizeof (nphi));
os.write ((char *) &ntheta, (char) sizeof (ntheta));
 Fkt.binWrite (os);
}

void Funsurf::binRead (ifstream &is)
{
 type=FUNSURF;
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
 Fkt.binRead(is);
}
