# cmake version
set(headers
    vector.h 
    vector_global.h 
    bessel.h 
    cfresnel.h 
    formio.h 
    fresnel.h 
    gamma.h 
    getreso.h 
    gitter.h 
    goodies.h 
    hoch.h 
    intersect.h
    matrix.h 
    resutil.h
    spherical.h 
    t.h 
    tunnel.h 
    vektor.h)

  set(vector_sources
    vector.cpp 
    bessel.cc 
    cfresnel.cc 
    fresnel.cc 
    gamma.cc 
    getreso.cc 
    gitter.cc 
    hoch.cc 
    io.cc 
    intersect.cpp
    tunnel.cc 
    goodies.cpp 
    matrix.cpp 
    spherical.cpp)

ADD_LIBRARY(vector ${vector_sources})
