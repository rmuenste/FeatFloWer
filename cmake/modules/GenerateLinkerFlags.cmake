# In this cmake file we configure
# the libraries and other parameters of the
# applications related to the linking process
#
#if(USE_OPTICALTWEEZERS)
#  add_subdirectory(extern/libraries/opticaltweezers)
#endif(USE_OPTICALTWEEZERS)

if(WIN32)
set(FF_DEFAULT_LIBS
  amd
  umfpack4
  feat2d
  feat3d
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  )
else(WIN32)
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
endif(WIN32)

if(WIN32)
  set(FF_APPLICATION_LIBS
    amd
    umfpack4
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    inshape3dcore
    Utility
    Math
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    ff_cinterface
    ff_postprocessing
    ff_ini_aux
    ff_ini_c
    ff_util
    ff_assemblies
    ff_mesh
    ff_elements
    ff_le_solvers
    ff_fbm
    ff_LinSc
    ff_q2p1
    ff_quadLS_app
    )
else(WIN32)
#  set(FF_APPLICATION_LIBS
#    amd
#    umfpack4
#    feat2d
#    feat3d
#    ${BLAS_LIBRARIES}
#    ${LAPACK_LIBRARIES}
#    stdc++
#    stdc++fs
#    inshape3dcore
#    Utility
#    Math
#    ${LIBRT_LIBRARY}
#    ${MPI_Fortran_LIBRARIES}
#    )
  set(FF_APPLICATION_LIBS
    amd
    umfpack4
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    stdc++
    stdc++fs
    inshape3dcore
    Utility
    Math
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    ff_cinterface
    ff_postprocessing
    ff_ini_aux
    ff_ini_c
    ff_util
    ff_assemblies
    ff_mesh
    ff_elements
    ff_le_solvers
    ff_fbm
    ff_LinSc
    ff_q2p1
    ff_quadLS_app
    )
endif(WIN32)


#${Elements} 
#${src_q2p1} 
#${src_pp3d} 
#${src_mpi} 
#${src_PLin} 
#${src_LinSc} 
#${src_quadLS_app} 
#${src_visco} 
#${src_mesh}
#${src_cinterface}
#${postprocessing}


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

if(USE_BOOST)
  list(APPEND FF_APPLICATION_LIBS ${Boost_LIBRARIES})
endif(USE_BOOST)
  
endif(USE_CGAL)

if(USE_OPTICALTWEEZERS)

  list(APPEND FF_APPLICATION_LIBS ${OPTICALTWEEZERS_LIBRARIES})
  
endif(USE_OPTICALTWEEZERS)


if(USE_HYPRE)

  list(APPEND FF_DEFAULT_LIBS ${HYPRE_LIBRARIES})

  list(APPEND FF_APPLICATION_LIBS ${HYPRE_LIBRARIES})

endif(USE_HYPRE)
 


