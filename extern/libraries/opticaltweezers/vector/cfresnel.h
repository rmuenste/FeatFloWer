#ifndef dc
#define dc double_complex
#endif

#ifndef FRESNEL_H
#define FRESNEL_H
#include <complex>

#ifndef double_complex
#define double_complex complex<double>
#endif double_complex

#include "vector.h"
#include "resutil.cc"


double_complex Fresnel_trans (double_complex alpha, double n1, double n2);
#endif
