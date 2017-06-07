#include <iostream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include "hoch.h"
#include "t.h"
#include "gamma.h"
void jnyp (double x, double ny, int g, double * jny)
/****************************************************************************
C**** JNYP BERECHNET DIE SPHAERISCHEN ZYLINDERFUNKTIONEN MITTELS RUECK-        *
C**** WAERTSREKURSION (SIEHE GAUTSCHI, SIAM-REVIEW, VOL. 9, NO. 1, JAN 1967)   *
C**** DER 4. PARAMETER ENTHAELT NACH AUFRUF DER ROUTINE DIE WERTE VON          *
C**** JNY-1(X) UND JNY(X).                                                     *
C*******************************************************************************/
{
double  j1[3], j2[3], rn[201], jnyh[3], nya, a; 
double eps, su, d, h1, h2, r, s, la, lah, f;  
int  m, n, i, j, ifail;
 
bool abfra;
j1[0] = 0.0;
j1[1] = 0.0;
j1[2] = 0.0;
ifail = 0;

nya = ny + 1.5;
m = (int)nya;
a = nya - m;
eps = 0.5 * hoch (10 , -g);
su = exp (a * log (0.5 * x));
d = 2.3026 * g + 1.3863;
h1 = 1.359 * x + t(0.7358 * ( d / x));
h2 = m * t(0.5 * d * m);
    if (h1 >  h2)
        {
         n = 1 + (int)h1;
	 }

    else
        {
         n=1 + (int)h2;  
         }   

abfra = true;

j = 1;
while ((j <= 10) && (abfra || (j <=  2)))
{
 lah = 1;

 for (i = 1; i <= n / 2; i = i++)
 {
 lah = lah * (a + i) / (i + 1);
 }
 r = 0;
 s = 0;

 for (i = n; i >= 1; i = i-1)
 {
 r = x / (2.0 * (a + i) -x * r);
 
	if ((i % 2) == 1) 
	{
	la = 0;
	}
	else 
	{
	lah = lah * (i + 2.0) / (i + 2.0 * a);
	la = (a + i) * lah;
	}
 
	s = r * (la + s);
   if (i <=  m ) 
   {
       rn[i] = r;
    }}
 f = su / (1 + s);
    for (i = 1; i<= m - 1; i = i++)
 {
 f = f * rn [i];
 }
 j2[1] = f;
 j2[2] = rn[m] * f;
 abfra = false;
 i = 1;
 while  ((! abfra) && (i <= 2)) 
     {
     abfra = fabs(j1[i] - j2[i]) > (eps * fabs (j2[i]));
     i = i + 1 ;
     }
 j1[1] = j2[1];
 j1[2] = j2[2];
 n = n + 5;
 j = j + 1;

}

jnyh[2] = j2[1];
jnyh[1] = (2 * ny + 1.0) / x * j2[1] - j2[2];
j2[1] = jnyh[1];
j2[2] = jnyh[2];
 r = sqrt (M_PI / (2.0 * x)) / Gamma(1.0 + a);
 /*cout << "r=" << r << endl;*/
jny[1] = j2[1] * r;
jny[2] = j2[2] * r;
}



/****************************************************************************
C*****  T IST EINE HILFSFUNKTION ZUR BERECHNUNG DES ANFANGS DER RUECK-         *
C*****  WAERTSREKURSION                                                        *
C*******************************************************************************/

double t (double x)
{
 double  l ;
 double t0, t1, t2, t3, t4, t5;
 double Erg;

 t0 = 1.0125;
 t1 = 0.8577;
 t2 = -0.129013;
 t3 = 0.0208645;
 t4 = -0.00176148;
 t5 = 0.000057941;

if (x <= 10)
 {
 Erg = t0 + x * (t1 + x * (t2 + x * (t3 + x * (t4 + x * t5))));
 }
else
{
l = log(x) - 0.775;
Erg = (0.775 - log(l)) / (1 + l);
Erg = 1.0 / (1.0 + Erg);
Erg = x / (l * Erg);
}
 return Erg;
}




/******************************************************************************
C****  JNYN LIEFERT JNY(X) UND JNY+1(X) ZU NEGATIVEM INDEX NY                 *
C*****************************************************************************/

void jnyn (double x, double ny, int g, double *jny)
{
 double j1[3], nya;
 int m, i;

 m = 1 + (int) (-ny);
 nya = ny + m;
 jnyp(x, nya + 1.0, g, jny);
 for (i = 1; i <= m; i = i++ )
 {
   j1[1] = jny [1];
   j1[2] = jny [2];
   jny[1] = (2.0 * nya + 1.0) * j1[1] / x - j1[2];
   jny[2] = j1[1];
   nya = nya - 1;
  }
}


/*
   besseldiff liefert die 1. Ableitung der sph. Besselfkt nach Rekursion
*/

void besseldiff (double x, double ny, int g, double &jny, double &jnydiff)

{
 double *jn;

 jn = new double[3];
 jnyp(x,ny,g,jn);
 jny  = jn[2];
 jnydiff = jn[1] - (ny+1)/x*jn[2];
delete jn;
}

/*
   JNYN LIEFERT JNY(X) UND JNY+1(X) ZU NEGATIVEM INDEX NY
*/

void jnynn (double x, double ny, int g, double * jny)
{
 double jnn[3], j1[3], nya;
 int  m, i;

 m = 1 + (int)-ny;
 nya = ny + m;

 jnyp (x, nya + 1, g, jnn);

 for (i = 1; i <= 800; i = i++)
 {
     j1[1] = jnn [1];
     j1[2] = jnn[2];
     jnn[1]= (2 * nya + 1) * j1[1] / x - j1[2];
     jnn[2] = j1[1];
     nya = nya - 1;
 }
}


/***  HNY2DIF LIEFERT HNY2(X) UND 1/X * D/DX (X * HNY2) **/
/*
  jetzt liefert es wirklich die Ableitung !!!
 */

void hny2dif (double x, double ny, int g, complex<double>  *hny1,
              complex<double>  *hny2)
{
 complex<double>  h[3], f1, f2;
 double jp[3], jn[3];
    
 
 f1 = complex<double>  (0.0, -1.0 / cos (ny * M_PI));
 f2 = complex<double>  (cos((ny + 0.5) * M_PI), sin((ny + 0.5) * M_PI));

 jnyp(x, ny, g, jp);
 jnyn (x, -ny - 1, g, jn);

 h[1] = f1 * (f2 * jp[2] - jn [1]);
 h[2] = f1 * (f2 * jp[1] + jn [2]);
 hny2[1] = h[1];
 hny2[2] = h[2] - (ny / x) * h[1];
 //hny1[2] = hny1[2] - hny1[1] / x;
 hny2[2] = hny2[2] - hny2[1] / x;
 hny1[1] = conj(hny2[1]);
 hny1[2] = conj(hny2[2]);
}
