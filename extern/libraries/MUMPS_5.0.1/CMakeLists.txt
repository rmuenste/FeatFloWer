project(MUMPS C CXX Fortran)

cmake_minimum_required(VERSION 3.1)

#add_executable(eigen-test eigen-test.cpp)

#set(CMAKE_MODULE_PATH ${CMAKE_ROOT}/Modules "${CMAKE_SOURCE_DIR}/cmake")
#message(STATUS "Module path ${CMAKE_MODULE_PATH}")

#target_include_directories(eigen-test PUBLIC ${CMAKE_SOURCE_DIR}/../eigen3)
#target_compile_options(eigen-test PUBLIC -std=c++11)
add_subdirectory(PORD)

#set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS})
#message(STATUS "CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS}")

set(mumps_common_src
src/ana_omp_m.F
src/tools_common_mod.F
src/mumps_static_mapping.F
src/mumps_sol_es.F
src/fac_future_niv2_mod.F
src/mumps_comm_ibcast.F
src/mumps_ooc_common.F
src/double_linked_list.F
src/fac_asm_build_sort_index_m.F
src/fac_asm_build_sort_index_ELT_m.F
src/omp_tps_common_m.F
src/mumps_l0_omp_m.F
src/front_data_mgt_m.F
src/fac_maprow_data_m.F
src/fac_descband_data_m.F
src/fac_ibct_data_m.F
src/ana_orderings.F
src/ana_AMDMF.F
src/bcast_errors.F
src/estim_flops.F
src/mumps_type_size.F
src/mumps_type2_blocking.F
src/mumps_version.F
src/tools_common.F
src/mumps_print_defined.F
src/mumps_common.c
src/mumps_orderings.c
src/mumps_size.c
src/mumps_io.c
src/mumps_io_basic.c
src/mumps_io_thread.c
src/mumps_io_err.c
src/ana_set_ordering.F
src/mumps_numa.c
)

add_library(mumps_common ${mumps_common_src})

set_target_properties(mumps_common PROPERTIES LINKER_LANGUAGE Fortran)

target_include_directories(mumps_common PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/PORD/include /sfw/openmpi/intel16.0.x/1.10.2/nonthreaded/include)

target_compile_definitions(mumps_common PRIVATE Add_=1 pord=1 MUMPS_ARITH=MUMPS_Arith_d)

target_link_libraries(mumps_common pord) 

set(dmumps_src
src/dana_aux_par.F
src/dana_lr.F
src/dfac_asm_master_m.F
src/dfac_asm_master_ELT_m.F
src/domp_tps_m.F
src/dmumps_comm_buffer.F
src/dmumps_load.F
src/dmumps_lr_data_m.F
src/dmumps_ooc_buffer.F
src/dmumps_ooc.F
src/dmumps_struc_def.F
src/dstatic_ptr_m.F
src/dlr_core.F
src/dfac_lr.F
src/dfac_omp_m.F
src/dfac_par_m.F
src/dfac_front_LU_type1.F
src/dfac_front_LU_type2.F
src/dfac_front_LDLT_type1.F
src/dfac_front_LDLT_type2.F
src/dfac_front_aux.F
src/dfac_front_type2_aux.F
src/dini_driver.F
src/dana_driver.F
src/dfac_driver.F
src/dsol_driver.F
src/dend_driver.F
src/dana_aux_ELT.F
src/dana_aux.F
src/dana_dist_m.F
src/dana_LDLT_preprocess.F
src/dana_reordertree.F
src/darrowheads.F
src/dbcast_int.F
src/dfac_asm_ELT.F
src/dfac_asm.F
src/dfac_b.F
src/dfac_distrib_distentry.F
src/dfac_distrib_ELT.F
src/dfac_lastrtnelind.F
src/dfac_mem_alloc_cb.F
src/dfac_mem_compress_cb.F
src/dfac_mem_free_block_cb.F
src/dfac_mem_stack_aux.F
src/dfac_mem_stack.F
src/dfac_process_band.F
src/dfac_process_blfac_slave.F
src/dfac_process_blocfacto_LDLT.F
src/dfac_process_blocfacto.F
src/dfac_process_bf.F
src/dfac_process_end_facto_slave.F
src/dfac_process_contrib_type1.F
src/dfac_process_contrib_type2.F
src/dfac_process_contrib_type3.F
src/dfac_process_maprow.F
src/dfac_process_master2.F
src/dfac_process_message.F
src/dfac_process_root2slave.F
src/dfac_process_root2son.F
src/dfac_process_rtnelind.F
src/dfac_root_parallel.F
src/dfac_scalings.F
src/dfac_determinant.F
src/dfac_scalings_simScaleAbs.F
src/dfac_scalings_simScale_util.F
src/dfac_sol_pool.F
src/dfac_type3_symmetrize.F
src/dini_defaults.F
src/mumps_c.c
src/dmumps_driver.F
src/dmumps_f77.F
src/dmumps_iXamax.F
src/dana_mtrans.F
src/dooc_panel_piv.F
src/drank_revealing.F
src/dsol_aux.F
src/dsol_bwd_aux.F
src/dsol_bwd.F
src/dsol_c.F
src/dsol_fwd_aux.F
src/dsol_fwd.F
src/dsol_matvec.F
src/dsol_root_parallel.F
src/dtools.F
src/dtype3_root.F
src/dmumps_restart.F)

add_library(dmumps ${dmumps_src})

set_target_properties(dmumps PROPERTIES LINKER_LANGUAGE Fortran)

target_include_directories(dmumps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${MPI_Fortran_INCLUDE_PATH})

target_compile_definitions(dmumps PRIVATE Add_=1 pord=1 MUMPS_ARITH=MUMPS_Arith_d)

target_link_libraries(dmumps mumps_common pord) 

