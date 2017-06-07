#include "funktion.h"
#include <strstream>
#include <string.h>
#ifndef min
#define min(a,b) a<b ? a : b
#endif

#ifndef max
#define max(a,b) a>b ? a : b
#endif

bool findeZahl (string s,int pos, int &lpos, int &rpos)
{
 lpos=s.find_last_not_of ("0123456789."); 
 rpos=s.find_last_not_of ("0123456789.e");
 // noch nicht fertig ! 
 return true; 
}

string boolout (bool x)
{
 return x ? string("true") : string("false");  
} 

// using namespace std;
Funktion :: Funktion(string fkt)
{
this->fktm=fkt;
this->S[0]=fkt;
 nS=1;
//pass(fktm,this->fkt,S,nS);
 pass();
} 

double  Funktion::operator ()(double x,const double y, const double z)
{
 return  interpret(x,y,z,S[0],S,nS);
}

void Funktion::pass ()
{
int i = 1;
signed int kz;
signed int ka ;
char ch[2];
char Str[20];
string str; 
signed int c,l,pos,pose;
S[0]=fktm;
nS=1;
l= S[0].length();
c=0;
ch[1]='\0'; 
while (c<l) // Ersetze Funktion durch Token
{
 for (int i=0; i<NFUNC; i++)
 {
  pos=S[0].find(func[i],c);
    ch[0]=token[i];
    if (pos>=0) S[0].replace(pos,strlen(func[i]),ch);
 }
 c++;
}
do
{
 pose=S[0].find_first_of ("eE");
 if (pose>0)
 {
  S[0].insert(pose+1,"(");
  S[0][pose]='k';
  str=S[0].substr(pose+3,S[0].length()-pose);
  pos=str.find_first_not_of ("(1234567890.k");
  if (pos==-1) pos=S[0].length();
  else pos=pose+3+pos;

   S[0].insert(pos,")");
  S[0].insert(pose,"*");
 }
} while (pose>0);


do
{
 kz=S[0].find_first_of(")");
 if (kz>=0)
 {
  ka=S[0].rfind("(",kz);
  S[i]=S[0].substr(ka+1,kz-ka-1);
  sprintf (Str,"S%i",i);
  S[0]=S[0].replace(ka,kz-ka+1,Str);
  i=i+1;
 }
}
while(kz>0);
nS=i;
}


double Funktion::interpret (double x, double y, double z, string f, string *S, int nS)
{
 double erg;
 signed int pos;
 string h;
 pos=f.find_first_of("+");
 if (pos==0) return interpret(x,y,z,f.substr(pos+1),S,nS); 
 else
 if (pos>0) return interpret(x,y,z,f.substr(0,pos),S,nS)+interpret(x,y,z,f.substr(pos+1),S,nS); 
 pos=f.find_last_of("-");
 if (pos==0) return -interpret(x,y,z,f.substr(1),S,nS);
 if (pos>0) return interpret(x,y,z,f.substr(0,pos),S,nS)-interpret(x,y,z,f.substr(pos+1),S,nS);
 pos=f.find_first_of("*");
 if (pos>=0) return interpret(x,y,z,f.substr(0,pos),S,nS)*interpret(x,y,z,f.substr(pos+1),S,nS);
 pos=f.find_first_of("/");
  if (pos>=0) return interpret(x,y,z,f.substr(0,pos),S,nS)/interpret(x,y,z,f.substr(pos+1),S,nS);
 pos=f.find_first_of("^");
 if (pos>=0)
  {
   return pow (interpret(x,y,z,f.substr(0,pos),S,nS),interpret(x,y,z,f.substr(pos+1),S,nS));
  }
 switch (f[0])
 {
  case ZEHNHOCH : return pow(10.0,interpret (x,y,z,f.substr(1),S,nS)); break;
  case ASIN : return asin(interpret (x,y,z,f.substr(1),S,nS)); break;
  case ACOS : return acos(interpret (x,y,z,f.substr(1),S,nS)); break;
  case ATAN : return atan(interpret (x,y,z,f.substr(1),S,nS)); break;
  case SIN  : return sin(interpret (x,y,z,f.substr(1),S,nS)); break;
  case COS  : return cos(interpret (x,y,z,f.substr(1),S,nS)); break;
  case TAN  : return tan(interpret (x,y,z,f.substr(1),S,nS)); break;
  case LN   : return log(interpret (x,y,z,f.substr(1),S,nS)); break;
  case EXP  : return exp(interpret (x,y,z,f.substr(1),S,nS)); break;
  case SQRT : erg= sqrt(interpret (x,y,z,f.substr(1),S,nS)); return erg; break;
  case 'x' : return x; break;
  case 'y' : return y; break;
  case 'z' : return z; break;
  case 'S' : char *Str;
             h= f.substr(1);
             return interpret (x,y,z,S[atoi(h.c_str())],S,nS); break;
  case '0' :
  case '1' :
  case '2' :
  case '3' :
  case '4' :
  case '5' :
  case '6' :
  case '7' :
  case '8' :
  case '9' : return atof(f.c_str()); break;
 }


}

string diff (string f, Funktion &F)
{
   char hStr[200], hStr2[200];
   string h1,q;
   signed int pos,t,n;
   int l,ns;
 string h,f1,f2;
 pos=f.find_first_of("+");
  if (pos>=0) {
               f1=diff (f.substr(0,pos),F);
               f2=diff (f.substr(pos+1),F);
                              h=f1+"+"+f2;
                  return h;}
 pos=f.find_last_of("-");
 if ((pos>0) && (f[pos-1]!='^')) {
             if (pos>=f.length()-1) l=f.length()-2;
             else l=pos+1; 
              h=diff (f.substr(0,pos),F)+"-"+diff (f.substr(l),F);
             return h;}
 else if (pos==0)
            { 
             ns=F.nS;
             if (F.nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             F.nS++;
             F.S[ns]="-"+diff (f.substr(pos+1),F);
             F.clean(); 
             sprintf (hStr,"S%i",ns);
             h=hStr;
             return h;
            }

 pos=f.find_first_of("*");
 if (pos>=0) {
             f1=f.substr(0,pos);
             f2=f.substr(pos+1,f.length()-pos);
             ns=F.nS; 
             sprintf (hStr,"S%i*",F.nS);
             sprintf (hStr2,"+S%i*",F.nS+1);
             h="";
             h=string(hStr)+f2+string(hStr2)+ f1;
             F.nS+=2;  
             F.S[ns]=diff (f1,F);
             ns++;
             F.S[ns]=diff(f2,F); 
             return h;}



 pos=f.find_first_of("/");
 if (pos>=0) {
             int ns2;
             string hb;
              ns=F.nS;
             f1=f.substr(0,pos);
             f2=f.substr(pos+1);
             if (F.nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             F.nS++;
             hb=diff (f1,F)+"*"+f2+"-"+f1+"*"+diff(f2,F);
             F.S[ns]=hb; //diff (f1,F)+"*"+f2+"-"+f1+"*"+diff(f2,F);
             if (F.nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             ns2=F.nS;
             F.nS++;
             F.S[ns2]=f2+"*"+f2;
             char hStr[200];
             sprintf (hStr,"S%i/S%i",ns,ns2);
             h=hStr;
             return h; }

 pos=f.find_first_of("^");
 if (pos>=0)
          {
            f1=f.substr(0,pos);
            f2=f.substr(pos+1);
            switch (f2[0])
            {
             case '-' :
             case '0' :
             case '1' :
             case '2' :
             case '3' :
             case '4' :
             case '5' :
             case '6' :
             case '7' :
             case '8' :
             case '9' : t=f2.find(".");
                        if (t<0) // ganze Zahl => Typ: f1^n n ganz
                        {
                         n=(signed int)atof (f2.c_str());
                         if (n-1>=0)  
                         {  
                          sprintf (hStr,"%i*%s^%i*",n,f1.c_str(),n-1);
                          ns=++F.nS;
                          F.S[ns]=diff(f1,F);
                          h=hStr;
                          sprintf (hStr,"S%i",ns);
                          h=h+hStr;
                          F.nS++;    
                         }
                         else 
                          {
                           sprintf (hStr,"%i",n-1); 
                           F.S[F.nS]=hStr;
                           sprintf (hStr,"%i*%s^S%i",n,f1.c_str(),F.nS);
                           h=hStr;
                           F.nS++; 
                          }  
                          return h;
                        }
             }
             //h=f+"*"+diff (f2,F)+"*"+LN+f1+"+"+f+"*"+f2+"*"+diff(f1,F)+"/"+f1;  // ?????????????
             h="überraschung";
             return h;

          }

switch (f[0])
 {
  case SIN : f1=f.substr(1);  // Ableitung von sin
             h=COS+f1+"*S";
             if(F.nS>=SMAX) return "ERROR";
             ns=F.nS++;
             F.S[ns]=diff(f1,F);
             sprintf(hStr,"%i",ns);
             h=h+hStr;
             return h; break;
 
  case EXP : f1=f.substr(1);  // Ableitung von sin
             h=EXP+f1+"*S";
             if(F.nS>=SMAX) return "ERROR";
             ns=F.nS++;
             F.S[ns]=diff(f1,F);
             sprintf(hStr,"%i",ns);
             h=h+hStr;
             return h; break;

 
   case COS : f1=f.substr(1);  // Ableitung von cos 
             h=string("-")+SIN+f1+"*S";
             if(F.nS>=SMAX) return "ERROR";
             ns=F.nS++;
             F.S[ns]=diff(f1,F);
             sprintf(hStr,"%i",ns);
             h=h+hStr;
             return h; break;
 
  case TAN : f1=f.substr(1);  // Ableitung von tan=1/cos^2 
             ns=F.nS++;
             sprintf(hStr,"1/S%i^2",ns);
             h=hStr;
             F.S[ns]=COS+f1;
             if(F.nS>=SMAX) return "ERROR";
             ns=F.nS++;
             F.S[ns]=diff(f1,F);
             sprintf(hStr,"%i",ns);
             h=h+"*S"+hStr;
             return h; break;

  case LN : f1=f.substr(1);  // Ableitung von ln=1/x 
             ns=F.nS++;
             sprintf(hStr,"1/S%i",ns);
             h=hStr;
             F.S[ns]=f1;
             if(F.nS>=SMAX) return "ERROR";
             ns=F.nS++;
             F.S[ns]=diff(f1,F);
             sprintf(hStr,"%i",ns);
             h=h+"*S"+hStr;
             return h; break;
 
  case SQRT :  
             f1=f.substr(1);
             if (F.nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             F.S[F.nS]="-0.5*1/i"+f1;
             F.nS++;
             sprintf(hStr,"S%i*",F.nS-1);
             h=hStr+diff (f1,F);
             return h; break;

  case 'x' :  h="1"; return h; break;
//  case 'x' : 
//  case 'y' : 
//  case 'z' : c=(int)(f[0]==var)+'0'; h=c;  return h; break;
 case 'S' : h=f.substr(1);
             h1=diff (F.S[atoi(h.data())],F); 
             return h1;
             break;
             // F.S[atoi(h.data())]=diff (F.S[atoi(h.data())],F); 
//             return F.S[atoi(h.data())]; 

//             return diff (S[atoi(h.data())],S,nS); break;
  case '#':
  case 'y':
  case 'z':
  case '0' :
  case '1' :
  case '2' :
  case '3' :
  case '4' :
  case '5' :
  case '6' :
  case '7' :
  case '8' :
  case '9' : h="0"; return h; break;

 }
}
/* string diff (string f,  string *S, int &nS)
{
    char var;
   char c;
   char hStr[200];
   signed int pos,t,n;
 string h,f1,f2;
 pos=f.find_first_of("+");
  if (pos>=0) {
               f1=diff (f.substr(0,pos),S,nS);
               f2=diff (f.substr(pos+1),S,nS);
                              h=f1+"+"+f2;
                  return h;}
 pos=f.find_last_of("-");
 if (pos>0) {
              h=diff (f.substr(0,pos),S,nS)+"-"+diff (f.substr(pos+1),S,nS);
             return h;}
 if (pos==0)
            {
             if (nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             S[nS]="-"+diff (f.substr(pos+1),S,nS);
             nS++;
             sprintf (hStr,"S%i",nS-1);
             h=hStr;
             return h;
            }

 pos=f.find_first_of("*");
 if (pos>=0) {
             f1=f.substr(0,pos);
             f2=f.substr(pos+1);
             h=diff (f1,S,nS)+"*"+f2+"+"+ f1+"*("+diff (f2,S,nS)+")";
             return h;}



 pos=f.find_first_of("/");
 if (pos>=0) {
             f1=f.substr(0,pos);
             f2=f.substr(pos+1);
             if (nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             S[nS]=diff (f1,S,nS)+"*"+f2+"-"+f1+"*"+diff(f2,S,nS);
             nS++;
             if (nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             S[nS]=f2+"*"+f2;
             nS++;
             char hStr[200];
             sprintf (hStr,"S%i/S%i",nS-2,nS-1);
             h=hStr;
             return h; }

 pos=f.find_first_of("^");
 if (pos>=0)
          {
            f1=f.substr(0,pos);
            f2=f.substr(pos+1);
            switch (f2[0])
            {
             case '-' :
             case '0' :
             case '1' :
             case '2' :
             case '3' :
             case '4' :
             case '5' :
             case '6' :
             case '7' :
             case '8' :
             case '9' : t=f2.find(".");
                        if (t<0) // ganze Zahl => Typ: f1^n n ganz
                        {
                         n=(signed int)atof (f2.c_str());
                         if (n-1>=0)   
                          sprintf (hStr,"%i*%s^%i*",n,f1.c_str(),n-1);
                         else 
                          {
                           sprintf (hStr,"%i",n-1); 
                           S[nS]=hStr;
                           sprintf (hStr,"%s^S",n,f1.c_str(),nS);
                           nS++; 
                          }  
                          h=hStr+diff(f1,S,nS);
                          return h;
                         
                        }
             }
        h=f+"*"+diff (f2,S,nS)+"*"+LN+f1+"+"+f+"*"+f2+"*"+diff(f1,S,nS)+"/"+f1;  // ??????????????
             return h;

          }

switch (f[0])
 {
  case SIN : f1=f.substr(1);  // Ableitung von sin
             h="e"+f1+"*S";
             if(nS>=SMAX) return "ERROR";
             S[nS]=diff(f1,S,nS);
             nS++;
             sprintf(hStr,"%i",nS-1);
             h=h+hStr;
             return h; break;
  case COS : f1=f.substr(1);  // Ableitung von cos
             if (nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             S[nS]="-d"+f1;
             nS++;
             sprintf(hStr,"S%i*",nS-1);
             h=hStr+diff (f1,S,nS);
             return h; break;
  case TAN : f1=f.substr(1); // Ableitung von tan=1/cos^2
             if (nS>=SMAX) return "ERROR"; // Fehler : zu viele dS-Stacks
             S[nS]="e"+f1+"*"+"e"+f1;
             nS++;
             sprintf(hStr,"1/S%i*",nS-1);
             h=hStr+diff (f1,S,nS);
             return h;
  case EXP : f1=f.substr(1); h=f+"*"+diff(f1,S,nS); return h; break;
  case 'x' :  h="1"; return h; break;
//  case 'x' : 
//  case 'y' : 
//  case 'z' : c=(int)(f[0]==var)+'0'; h=c;  return h; break;
  case 'S' : h=f.substr(1);
             S[atoi(h.data())]=diff (S[atoi(h.data())],S,nS); 
             return S[atoi(h.data())]; 
             break;
//             return diff (S[atoi(h.data())],S,nS); break;
  case '0' :
  case '1' :
  case '2' :
  case '3' :
  case '4' :
  case '5' :
  case '6' :
  case '7' :
  case '8' :
  case '9' : h="0"; return h; break;

 }
}
*/
Funktion Funktion::partdiff(string p)
{
int posx,posp,posq;
Funktion df;
string str=fktm;
if(p!="x")
{
posx=0;
 while (posx>=0)
  {
 posx=str.find_first_of("x");
 if(posx>=0) str=str.replace(posx,1,"#");
 posx=str.find_first_of("x");
 }

posp=0;
while (posp>=0)
 {
posp=str.find_first_of(p);
if(posp>=0) str=str.replace(posp,1,"x");
posp=str.find_first_of(p);
 }
Funktion h=Funktion(str);
df=*this;
df=diff(h);
string dstr=df.S[0];
posx=0;
while (posx>=0)
 {
posx=df.S[0].find_first_of("x");
if(posx>=0) df.S[0]=df.S[0].replace(posx,1,p);
posx=df.S[0].find_first_of("x");
 }

posq=0;
while (posq>=0)
 {
posq=df.S[0].find_first_of("#");
if(posq>=0) df.S[0]=df.S[0].replace(posq,1,"x");
posq=df.S[0].find_first_of("#");
 }
posx=0;
while (posx>=0)
 {
posx=df.fktm.find_first_of("x");
if(posx>=0) df.fktm=df.fktm.replace(posx,1,p);
posx=df.fktm.find_first_of("x");
 }

posq=0;
while (posq>=0)
 {
posq=df.fktm.find_first_of("#");
if(posq>=0) df.fktm=df.fktm.replace(posq,1,"x");
posq=df.fktm.find_first_of("#");
 }


for(int i=1 ; i<df.nS; i++)  //ersetzt in den stacks
  {
 posx=0;
 while (posx>=0)
 {
posx=df.S[i].find_first_of("x");
if(posx>=0) df.S[i]=df.S[i].replace(posx,1,p);
posx=df.S[i].find_first_of("x");
 }

posq=0;
while (posq>=0)
 {
posq=df.S[i].find_first_of("#");
if(posq>=0) df.S[i]=df.S[i].replace(posq,1,"x");
posq=df.S[i].find_first_of("#");
 }


}//end for Stack

return df;
}


else
{
//Funktion f(fkt);
df=*this;
df=diff(df); 
}
//df.clean();
return df;
}

ostream & operator << (ostream & os,Funktion &f)
{ 
 os << "fktm=" << f.fktm << endl;
 os << "nS=" << f.nS << endl;
 for (int i=0; i<f.nS; i++)
 os << "S[" << i << "]=" << f.S[i] << endl;
 return os; 
}

void Funktion::print ()
{
 for (int i=0; i<nS; i++)
 {
  cout << "S" << i << "=";
  printStack(i);
 }  
}


void Funktion::printStack(int i)
{
 int n=S[i].length();
 for (int l=0; l<n; l++)
 {
  switch (S[i][l])
  { 
   case ASIN : cout << "asin"; break;
   case ACOS : cout << "acos"; break;
   case ATAN : cout << "atan"; break;
   case SIN : cout << "sin"; break;
   case COS : cout << "cos"; break;
   case TAN : cout << "tan"; break;
   case LN  : cout << "ln"; break;
   case EXP : cout << "exp"; break;
   case SQRT : cout << "sqrt"; break;
   case ZEHNHOCH : cout << "10^"; break;
   default : cout << S[i][l]; 
  }
 }
 cout << endl;   
}

Funktion &  Funktion:: operator = (const Funktion & f)
{
if( this== &f) return *this;

nS=f.nS;
fktm=f.fktm;
for(int i=0; i< f.nS; i++)
 this->S[i]=f.S[i];
return *this;
}



Funktion diff (Funktion &f)
{
 Funktion df=f;
// df.fkt=diff (df.fkt,df.S,df.nS); 
   df.S[0]=diff (df.S[0],df); 
  df.fktm=df.S[0];
//  df.pass();
//  df.clean();
 return df;
}


string Funktion::cleanStack(int j)
{
 //erst mal nach einer einzelnen Null suchen
 int h,rpos,lpos,pos=0;
 int lmpos,rmpos,hlpos,hrpos; 
 bool found=true,fnd;
 bool hoch=false;
 bool ok,stop;
 string s;
 s=S[j];
 pos=s.find_first_of('0');
 if (s.length()>1)  
 do 
 {
   found=pos>=0;
   if(found) 
   {
    for (int i=0; (i<NZAHL) && found; i++)  // erst mal schauen, ob die Null auch alleine steht 
      {
       found=(s[pos+1]!=zahl[i]);   // danach steht keine Ziffer
       found=found && ((s[pos-1]!=zahl[i]) || (pos==0));   // davor auch nicht
      }
 
    if (found) // da steht wirklich eine einsame 0
    {
     ok=true;
     hoch=false;
     
     lpos=s.find_last_of ("+-",pos);
     lmpos=s.find_last_of ("*",pos);
     rpos=s.find_first_of ("+-",pos);
     rmpos=s.find_first_of ("*",pos); 
     hlpos=max (lpos,lmpos);
     if ((rmpos<rpos) && (rmpos>=0)) hrpos=rmpos; else hrpos=rpos;
     if (hlpos<0) { lpos=0; }
     if (pos>0) { if (s[pos-1]=='^') s.replace (hlpos+1,pos-hlpos,"1"); else pos=-1;} // Da steht ein ^ vor der 0
     if (pos<=0) 
     {
      int ds; 
      if (rpos<0) rpos=s.length();
            if (lpos<=0) { lpos=0; ds=rpos+1;} else ds=rpos-lpos;
//      cout << "lpos=" << lpos << "   rpos=" << rpos << "   ds=" << ds << endl;

      if (s[rpos]=='-') s.erase (lpos,rpos-lpos);  // ein "-" steht rechts, das muss stehen bleiben
      else { s.erase (lpos,ds);}
     }
    } 
    //pos=-1;
   }
 pos=s.find('0',pos+1);
 }  while ((s.length()>1) && (pos>=0));

// jetzt mal nach den 1-sen suchen

 pos=0;
 h=s.length();
stop=false;
while ((pos>=0) && (h>=2) && !stop)
 {  
   stop=false; 
   pos=s.find_first_of("1",pos);
   found=(pos>=0) && (h>1);
   if(found) 
   {
    fnd=true; 
    for (int i=0; (i<NZAHL) && fnd; i++)  // erst mal schauen, ob die Eins auch alleine steht 
      {
       fnd=((s[pos+1]!=zahl[i]) ||(pos==s.length()-1 ));   // danach steht keine Ziffer
       fnd=fnd && ((s[pos-1]!=zahl[i]) || (pos==0));   // davor auch nicht
      }
     fnd = fnd && (s[pos-1]!='-'); // ist das vielleicht aber eine -1 ?
     fnd = fnd && (s[pos-1]!='S');  // oder ist das etwa Stack 1 ?
      if (fnd) // da steht wirklich eine einsame 1
      {
       if ((s[pos-1]=='^') || (s[pos-1]=='*')) { s.erase (pos-1,2);}
       else if ((s[pos+1]=='^') || (s[pos+1]=='*')){ s.erase (pos,2);} 
      }
      pos++;
   }
   else pos=-1;
  h=s.length();
  stop=pos>=h;
  }
 return s; 
}

void Funktion::clean ()
{
 bool found;
 int l,pos, i=0;
 string s; 
 char hstr[255];
 do  // Erst mal nach den Nullen suchen
 {  
    if ((S[i][0]=='0') && (S[i].length()==1))   // Der ganze Stack ist 0, muss also eliminiert werden
    {
     pos=-1;  
     sprintf (hstr,"S%i",i);
     int l=0;
     do
     {
       pos=S[l].find(hstr); 
       if (pos>=0) { S[l].replace(pos,strlen(hstr),"0"); l=-1; } // Link wurde gefunden und durch 0 ersetzt
       l++;
      }
      while (l<nS); 
 //    S[i]=""; 
      i=-1;
    }
    else 
    {
     S[i]=cleanStack(i);   // Stack i untersuchen

     if ((S[i][0]=='0') && (S[i].length()==1)) i=-1;
    }
  i++;
 }
 while (i<nS) ;
 
 if (S[0].length()==0) S[0]='0';
 i=0; 
 do  // Jetzt nach Einsen suchen
 {  
    S[i]=cleanStack(i);

    if ((S[i][0]=='1') && (S[i].length()==1) && (i>0))   // Der ganze Stack ist 1, muss also eliminiert werden
    {
     pos=-1;  
     sprintf (hstr,"S%i",i);
     l=0;
     do
     {
       pos=S[l].find(hstr);
       if (pos>=0) { S[l].replace(pos,strlen(hstr),"1"); l--; } // Link wurde gefunden und durch 1 ersetzt
       l++;
      }
      while (l<nS); 
     S[i]=""; 
     i=-1;
    }
    else 
    { 
     S[i]=cleanStack(i);   // Stack i untersuchen
     if ((S[i][0]=='1') && (S[i].length()==1) && (i>0)) i=-1;
    }
  i++;
 }
 while (i<nS) ;
 S[0]=cleanStack(0);
// if (S[0].length()==0) S[0]="0";
  
 // Jetzt schauen wir mal ob man einige Stacks direkt einsetzen können !
/* for (int i=1; i<nS; i++)
 {
  sprintf (hstr,"S%i",i);
  if (S[i].find_first_of ("+-")<0) 
  for (int j=0; j<nS; j++)
  {
   if ((S[j].length()!=0) && (i!=j))
   {
   pos=0;
   do
   { 
    found=false;
    pos=S[j].find(hstr,pos);
    if (pos>0) 
    {  
      for (int l=0; ((l<NFUNC) && (!found)); l++)  found=found || (S[j][pos-1]==token[l]);
      found=found || (S[j][pos-1]=='/'); //  || (S[j][pos-1]=='*');
    }
    if ((!found) && (pos>=0)) S[j].replace(pos,strlen(hstr),S[i]);
    pos++; 
   } while (pos>0);
  }
  } 
 }*/  
}

double Funktion::Nullstelle(double xstart, double xstop, double dx, bool &found, const double eps)
/* 
   Berechnet die nächste Nullstelle nach xstart bis xstop (Nur für Funktionen die
   nur von x abhängen (y,z=0)! 
   dx: grobe Suchweite (<Abstand zwischen 2 Nullstellen !)
   eps : Genauigkeit
*/
{
 double f,falt,fm;
 double x=xstart,xalt,xm;
 f=interpret(x,0,0,S[0],S,nS);
 // Erst mal grob die Nullstelle suchen 
 do
 {
  falt=f;
  xalt=x;
  x=x+dx;
  f=interpret(x,0,0,S[0],S,nS);
  found=(f*falt<0); 
 } while (!found && (x<xstop));
 if (!found) {/* cout << "nicht gefunden" << endl;*/ return xstart; }// Die Nullstelle wurde nicht gefunden

 // Jetzt suchen wir mal etwas genauer
 xm=x;
 do
 {
  xm=(x+xalt)/2.0;
  fm=interpret(xm,0,0,S[0],S,nS);
  if (fm*falt<0) { x=xm; f=fm; }
  else {xalt=xm; falt=fm;} 
 } while (x-xalt>eps);
 return xm;  
} 


void Funktion::binWrite (ofstream &os)
{
 os.write((char *) &fktm,(char) sizeof (fktm));
 os.write((char *) &nS, (char) nS);
 os.write((char *) &fktm,(char) sizeof (fktm));
 for (int i=0; i<nS; i++)
   os.write((char *) &S[i],sizeof(S[i]));
}

void Funktion::binRead (ifstream &is)
{
 is.read ((char *) &fktm,(char) sizeof (fktm));
 is.read ((char *) &nS, (char) nS);
 is.read ((char *) &fktm,(char) sizeof (fktm));
 for (int i=0; i<nS; i++)
   is.read ((char *) &S[i],sizeof(S[i]));
}


