# In this cmake file we configure
# the libraries and other parameters of the
# applications related to the linking process
#
#if(USE_OPTICALTWEEZERS)
#  add_subdirectory(extern/libraries/opticaltweezers)
#endif(USE_OPTICALTWEEZERS)

if(WIN32)
set(FF_DEFAULT_LIBS
  ${FF_UMFPACK_LINK_LIBS}
  feat2d
  feat3d
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${FF_MKL_PARDISO_LIBS}
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  )
else(WIN32)
set(FF_DEFAULT_LIBS
  ${FF_UMFPACK_LINK_LIBS}
  feat2d
  feat3d
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${FF_MKL_PARDISO_LIBS}
  ${LIBRT_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  stdc++fs
  )
endif(WIN32)

if(WIN32)
  set(FF_APPLICATION_LIBS
    ${FF_UMFPACK_LINK_LIBS}
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  ${FF_MKL_PARDISO_LIBS}
    inshape3dcore
    Utility
    Math
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    ff_cinterface
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
    ff_postprocessing
    )
else(WIN32)
  set(FF_APPLICATION_LIBS
    ${FF_UMFPACK_LINK_LIBS}
    feat2d
    feat3d
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  ${FF_MKL_PARDISO_LIBS}
    stdc++
    stdc++fs
    OpenMP::OpenMP_Fortran
    inshape3dcore
    Utility
    Math
    ${Boost_LIBRARIES}
    ${LIBRT_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
    ff_cinterface
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
    ff_postprocessing
    ff_particles
    )
endif(WIN32)

if(USE_PE)
  list(APPEND FF_DEFAULT_LIBS pe_static)
  list(APPEND FF_APPLICATION_LIBS pe_static)
endif(USE_PE)

if(USE_MUMPS)
  if(USE_EXTERNAL_MUMPS)
    set(MUMPS_LIBRARY_LIST
      ${FF_DMUMPS_LIBRARY} ${FF_MUMPS_COMMON_LIBRARY} ${FF_PORD_LIBRARY}
      ${FF_SCALAPACK_LIBRARY} ${FF_SCOTCH_LIBRARIES} ${FF_PARMETIS_LIBRARIES}
      ${FF_OPENBLAS_LIBRARY} gomp -lpthread -lm -ldl
      )
  else()
    set(MUMPS_LIBRARY_LIST
      dmumps mumps_common pord ${MKL_SCALAPACK_LIBRARY} -Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_BLACS_LIBRARY} -Wl,--end-group -lpthread -lm -ldl
      )
  endif()
  list(APPEND FF_DEFAULT_LIBS ${MUMPS_LIBRARY_LIST})
  list(APPEND FF_APPLICATION_LIBS ${MUMPS_LIBRARY_LIST})
endif(USE_MUMPS)

if(USE_ODE)
  list(APPEND FF_DEFAULT_LIBS ode)
  list(APPEND FF_APPLICATION_LIBS ode)
endif(USE_ODE)

if(USE_CGAL)
  list(APPEND FF_APPLICATION_LIBS "$<LINK_ONLY:CGAL::CGAL>")
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
 
