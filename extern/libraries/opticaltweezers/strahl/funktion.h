#ifndef FUNKTION_H
#define FUNKTION_H

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define SMAX 100

#define NFUNC 9 

#define ASIN 'a'
#define ACOS 'b'
#define ATAN 'c'
#define SIN  'd'
#define COS  'j'
#define TAN  'f'
#define LN   'g'
#define EXP  'h'
#define SQRT 'i'
#define ZEHNHOCH 'k'

#define NOPER 5 
#define PLUS '+'
#define MINUS '-'
#define MAL '*'
#define HOCH '^'

#define NZAHL 12 

const char func[NFUNC][5]={"asin","acos","atan","sin","cos","tan","ln","exp","sqrt"};

const char token[NFUNC]={ASIN,ACOS,ATAN,SIN,COS,TAN,LN,EXP,SQRT};
const char oper[NOPER]={'+','-','/','*','^'};
const char zahl[NZAHL]={'.','1','2','3','4','5','6','7','8','9','0','S'}; 
//const char oper[NOPER]={PLUS};
using namespace std;

class Funktion
{


// friend class Funsurf;
string cleanStack(int i);

void pass ();
void printStack(int i);
double interpret (double x, double y, double z, string f, string *S, int nS);

string fktm;
string S[SMAX];
int nS;
public:
double Nullstelle(double xstart, double xstop, double dx, bool &found, const double eps=1E-6);
void print();
void clean();
Funktion (){nS=0;};
Funktion (string fkt);
 //string getfkt() {return fkt; }
 string  getfktm() {return fktm; }
double  operator() (double x, const double y=0, const double z=0);
Funktion&  operator = (const Funktion & f);
Funktion partdiff(string p) ;
  void binWrite (ofstream &os);
  void binRead(ifstream &os); 

friend ostream& operator << (ostream &os, Funktion &f);
 friend Funktion diff (Funktion &f);
friend string diff (string f, Funktion &F);
};
#endif

   
