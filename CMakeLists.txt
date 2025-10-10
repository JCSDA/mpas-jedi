# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# MPAS
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( mpasjedi VERSION 3.1.0 LANGUAGES C CXX Fortran )

## Ecbuild integration
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )
option( ENABLE_MPAS_DOC "Build MPAS-JEDI documentation" OFF )

## Dependencies
find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS CXX Fortran )
find_package( Boost REQUIRED )

#OpenMP dependent packages
if(OpenMP_FOUND)
    if (MPAS_DOUBLE_PRECISION)
        find_package( MPAS 8.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere OpenMP)
    else()
        find_package( MPAS 8.0 REQUIRED COMPONENTS core_atmosphere OpenMP)
    endif()
    find_package( atlas 0.35.0 REQUIRED COMPONENTS OMP OMP_Fortran)
else()
    if (MPAS_DOUBLE_PRECISION)
        find_package( MPAS 8.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere)
    else()
        find_package( MPAS 8.0 REQUIRED COMPONENTS core_atmosphere)
    endif()
    find_package( atlas 0.35.0 REQUIRED )
endif()

# Use external ESMF if found
if(ESMF_FOUND)
  message(STATUS "Configure mpas-jedi for external ESMF")
  add_definitions(-DMPAS_EXTERNAL_ESMF_LIB)
  add_library(${PROJECT_NAME}::external::esmf ALIAS esmf)
else()
  message(STATUS "Configure mpas-jedi for internal ESMF")
endif()

find_package( rttov 12.1.0 QUIET )
find_package( ropp-ufo QUIET )
find_package( oops 1.10.0 REQUIRED )
find_package( saber 1.10.0 REQUIRED )
find_package( ioda 2.9.0 REQUIRED )
find_package( ufo 1.10.0 REQUIRED )

## Sources
add_subdirectory( src )

if( ENABLE_MPAS_DOC )
    add_subdirectory( docs )
endif()

## Tests
if( BUILD_TESTING )
  find_package( mpas-jedi-data QUIET )
  if( NOT mpas-jedi-data_FOUND )
    # Set values for downloading mpas-jedi test data from tarball
    set (mpas-jedi_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM/file")
    set (mpas-jedi_data_tag "3.1.0.jcsda")
    set (mpas-jedi_data_tarball_name "mpas-jedi_testinput_tier_1_${mpas-jedi_data_tag}.tar.gz")
    set (mpas-jedi_data_tarball_hash "1c9999dc04609f0c129ed39bc644bddc")
  endif()

  add_subdirectory( test )

  ecbuild_add_test( TARGET ${PROJECT_NAME}_coding_norms
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                    ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}/src
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
endif()

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
