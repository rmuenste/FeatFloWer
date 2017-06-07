#include "vector.h"
#include "fresnel.h"
#include "resutil.h"

complex<double>  tunnel (Vector<double> P, GlobalParms parms,
              double n1, double n2); 
//dc tunnel_reflect (const StrahlInfo &S, const GlobalParms &parms);
complex<double>  tunnel_reflect (GlobalParms parms, Vector<double> P, const StrahlInfo &S, double n1, double n2);



