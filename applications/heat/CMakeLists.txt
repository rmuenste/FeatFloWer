# CMakeLists of heat application

# source files for standad q2p1 
set(src_heat
app_init.f90
heat.f90
)

set(SSE_TARGET_NAME "heat")

# Add the application executable
add_executable(${SSE_TARGET_NAME} ${src_heat}) 

source_group(src_quadLS FILES ${src_q2p1})
source_group(src_pp3d FILES ${src_pp3d})
source_group(src_mpi FILES ${src_mpi})
source_group(src_PLin FILES ${src_PLin})
source_group(src_LinSc FILES ${src_LinSc})
source_group(src_quadLS FILES ${src_quadLS_ext})
source_group(src_visco FILES ${src_visco})
source_group(src_mesh FILES ${src_mesh})
source_group(src_cinterface FILES ${src_cinterface})
source_group(Elements FILES ${Elements})

if(GIT_FOUND)
  writeCommitHash(${SSE_TARGET_NAME})
endif(GIT_FOUND)

target_link_libraries(${SSE_TARGET_NAME} 
  ${FF_APPLICATION_LIBS}
  )

target_include_directories(${SSE_TARGET_NAME} PUBLIC ${FF_APPLICATION_INCLUDE_PATH})

target_compile_options(${SSE_TARGET_NAME} PRIVATE ${Fortran_FLAGS})

set_target_properties(${SSE_TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

add_dependencies(${SSE_TARGET_NAME} s3d_mesher)

get_target_property(MYTARGET ${SSE_TARGET_NAME} NAME)

IF(${OUT_OF_SOURCE_BUILD})
  include(${CMAKE_MODULE_PATH}/CreateDataDirectories.cmake)  
  createDefaultDirectories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF()

file(COPY heat_start.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY sensor_temperature_extraction.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
install(TARGETS ${MYTARGET}
    RUNTIME DESTINATION bin/${MYTARGET}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

set(DEF_DIRS start _data)
set(DEF_DIRS_BINARY _mesh _vtk _dump)

set(DEFAULT_CASE_FOLDER "")

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_adc/EWIKON_201912")
  list(APPEND DEFAULT_CASE_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/_adc/EWIKON_201912)
ENDIF()

FOREACH(dir ${DEF_DIRS})
  list(APPEND mydirs ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
ENDFOREACH()

FOREACH(dir ${DEF_DIRS_BINARY})
  list(APPEND mybindirs ${CMAKE_CURRENT_BINARY_DIR}/${dir})
ENDFOREACH()

install(DIRECTORY ${mydirs} DESTINATION bin/${MYTARGET})
install(DIRECTORY ${mybindirs} DESTINATION bin/${MYTARGET})

IF(DEFAULT_CASE_FOLDER)
  install(DIRECTORY ${DEFAULT_CASE_FOLDER} DESTINATION bin/${MYTARGET}/_ianus/HEAT)
ENDIF()


#=========================================================================
#                       INSTALL CONFIGURATION
#=========================================================================
# install the partitioner directory so the python sources can be found by
# the driver script
install(DIRECTORY ${CMAKE_SOURCE_DIR}/tools/partitioner DESTINATION bin/${MYTARGET})

set(HEAT_SCRIPTS heat_start.py)

list(APPEND HEAT_SCRIPTS sensor_temperature_extraction.sh)

install(FILES ${HEAT_SCRIPTS} DESTINATION bin/${MYTARGET}
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)

install(FILES "${CMAKE_SOURCE_DIR}/_data/MG.dat" DESTINATION bin/${MYTARGET}/_data 
    PERMISSIONS WORLD_EXECUTE OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE)
    