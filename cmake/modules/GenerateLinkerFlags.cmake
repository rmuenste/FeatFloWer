# In this cmake file we configure
# the libraries and other parameters of the
# applications related to the linking process
#
#if(USE_OPTICALTWEEZERS)
#  add_subdirectory(extern/libraries/opticaltweezers)
#endif(USE_OPTICALTWEEZERS)

set(FF_DEFAULT_LIBS
  amd
  umfpack4
  feat2d
  feat3d
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  stdc++fs
  )

if(WIN32)
  set(FF_APPLICATION_LIBS
    amd
    umfpack4
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    cdirs
    inshape3dcore
    Utility
    Math
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    )
else(WIN32)
  set(FF_APPLICATION_LIBS
    amd
    umfpack4
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    cdirs
    stdc++
    stdc++fs
    inshape3dcore
    Utility
    Math
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    )
endif(WIN32)

if(USE_MUMPS)

  set(MUMPS_LIBRARY_LIST
    dmumps mumps_common pord ${MKL_SCALAPACK_LIBRARY} -Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_BLACS_LIBRARY} -Wl,--end-group -lpthread -lm -ldl
    )

  list(APPEND FF_DEFAULT_LIBS ${MUMPS_LIBRARY_LIST})

  list(APPEND FF_APPLICATION_LIBS ${MUMPS_LIBRARY_LIST})

endif(USE_MUMPS)

if(USE_ODE)

  list(APPEND FF_DEFAULT_LIBS ode)

  list(APPEND FF_APPLICATION_LIBS ode)

endif(USE_ODE)

if(USE_CGAL)

  list(APPEND FF_APPLICATION_LIBS ${CGAL_LIBRARIES})

  list(APPEND FF_APPLICATION_LIBS ${Boost_LIBRARIES})
  
endif(USE_CGAL)

if(USE_OPTICALTWEEZERS)

  list(APPEND FF_APPLICATION_LIBS ${OPTICALTWEEZERS_LIBRARIES})
  
endif(USE_OPTICALTWEEZERS)

 


