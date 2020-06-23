include(findHelpers)
if(DEFINED ENV{RADAREMUL_VER})
  set(RADAREMUL_VER $ENV{RADAREMUL_VER})
  STRING(REGEX REPLACE " " "" RADAREMUL_VER ${RADAREMUL_VER})
endif()
if(NOT BUILD_RADAREMUL )
  if(DEFINED ENV{RADAREMUL_LIB})
    set(RADAREMUL_LIBRARY $ENV{RADAREMUL_LIB} )
  else()
    find_library( RADAREMUL_LIBRARY
      NAMES libradaremul.a 
      HINTS $ENV{COREPATH}/lib /usr/local/jcsda/nwprod_gdas_2014/lib
          ${COREPATH}/radaremul/v${RADAREMUL_VER}
          ${COREPATH}/radaremul/v${RADAREMUL_VER}/intel
          ${COREPATH}/radaremul/v${RADAREMUL_VER}/ips/${COMPILER_VERSION}
      PATH_SUFFIXES
        lib
       ${NO_DEFAULT_PATH}
      )
    message("Found RADAREMUL library ${RADAREMUL_LIBRARY}")
  endif()
endif()
if( NOT RADAREMUL_LIBRARY ) # didn't find the library, so build it from source
    message("Could not find RADAREMUL library, so building from libsrc")
    if( DEFINED ENV{RADAREMUL_SRC} )
      set( RADAREMUL_DIR $ENV{RADAREMUL_SRC} CACHE STRING "RADAREMUL Source Directory" )
    else()
      findSrc( "radaremul" RADAREMUL_VER RADAREMUL_DIR )
      set(RADAREMULINC  "${CMAKE_BINARY_DIR}/include")
    endif()
    set( radaremul "radaremul")
    set( BUILD_RADAREMUL "ON" CACHE INTERNAL "Build radaremul library" )
    #add_subdirectory(${CMAKE_SOURCE_DIR}/src/caps/radaremul)
    set( RADAREMUL_LIBRARY ${radaremul} )
    if( CORE_BUILT )
      list( APPEND CORE_BUILT ${RADAREMUL_LIBRARY} )
    else()
      set( CORE_BUILT ${RADAREMUL_LIBRARY} )
    endif()
else( NOT RADAREMUL_LIBRARY )
  if( CORE_LIBRARIES )
    list( APPEND CORE_LIBRARIES ${RADAREMUL_LIBRARY} )
  else()
    set( CORE_LIBRARIES ${RADAREMUL_LIBRARY} )
  endif()
endif( NOT RADAREMUL_LIBRARY )

set( RADAREMUL_LIBRARY_PATH ${RADAREMUL_LIBRARY} CACHE STRING "RADAREMUL Library Location" )
