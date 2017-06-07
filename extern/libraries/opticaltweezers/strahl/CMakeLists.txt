#-------------------------------------------------
#
# Project created by QtCreator 2016-05-19T13:06:21
#
#-------------------------------------------------
# output the system name

include_directories(../vector .)

set(strahl_sources
    strahl.cpp
    compound.cpp
    types.cpp
    zylinder.cpp
    ebene.cpp
    ellipsoid.cpp
    erythrocyte.cpp
    form.cpp
    funktion.cpp
    funsurf.cpp
    hohlfaser.cpp
    istrahl.cpp
    kegel.cpp
    kegelstumpf.cpp
    kegelstumpfhohl.cpp
    lgsurface.cpp
    LightSrc.cpp
    linse.cpp
    misc.cpp
    LightSrc.cpp
    linse.cpp
    misc.cpp
    PStrahl.cpp
    superellipsoid_n.cpp
    superellipsoid.cpp
    surface.cpp
    resutil.cpp
    koef.cc 
    mathfun.cpp 	
  )

set(headers_strahl
    strahl.h
    berechnung.h
    compound.h
    config.h    
    dreiecklist.h
    ebene.h
    ellipsoid.h
    erythrocyte.h
    form.h
    funktion.h
    funsurf.h
    hohlfaser.h
    istrahl.h
    kegel.h
    kegelstumpf.h
    kegelstumpfhohl.h
    lgsurface.h
    LightSrc.h
    linse.h
    misc.h
    PStrahl.h
    superellipsoid_n.h
    superellipsoid.h
    surface.h
    types.h
    zylinder.h
    resutil.h
    koef.h
    mathfun.h
  )

ADD_LIBRARY(strahl ${strahl_sources})

target_link_libraries(strahl vector)

