#ifndef __bessel__
#define __bessel__
#include <complex>

#ifndef double_complex
#define double_complex complex<double>
#endif double_complex

void jnyp  (double x, double ny, int g, double *jny);
void hny2dif (double x, double ny, int g, double_complex *hny1, 
              double_complex *hny2);
void besseldiff (double x, double ny, int g, double &jny, double &jnydiff);
#endif
