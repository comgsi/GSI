cmake_minimum_required(VERSION 2.6)

if(BUILD_RADAREMUL)
  file(GLOB RADAREMUL_SRC ${RADAREMUL_DIR}/*.f90 )
  set_source_files_properties( ${RADAREMUL_SRC} COMPILE_FLAGS ${RADAREMUL_Fortran_FLAGS} )
  add_library( ${radaremul} STATIC ${RADAREMUL_SRC} )
  set_target_properties( ${radaremul} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY} )

endif()
