# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# MPAS 
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( mpasjedi VERSION 1.0.0 LANGUAGES C CXX Fortran )

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
        find_package( MPAS 7.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere OpenMP)
    else()
        find_package( MPAS 7.0 REQUIRED COMPONENTS core_atmosphere OpenMP)
    endif()
    find_package( atlas 0.23.0 REQUIRED COMPONENTS OMP OMP_Fortran)
else()
    if (MPAS_DOUBLE_PRECISION)
        find_package( MPAS 7.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere)
    else()
        find_package( MPAS 7.0 REQUIRED COMPONENTS core_atmosphere)
    endif()
    find_package( atlas 0.23.0 REQUIRED )
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
find_package( oops 1.0.0 REQUIRED )
find_package( saber 1.0.0 REQUIRED )
find_package( ioda 1.0.0 REQUIRED )
find_package( ufo 1.0.0 REQUIRED )
# SACA
find_package( saca 1.0.0 QUIET )
if( ${saca_FOUND} )
    message(STATUS "SACA FOUND; Including SACA")
else( ${saca_FOUND} )
    message(STATUS "SACA NOT FOUND; Excluding SACA")
endif( ${saca_FOUND} )

## Sources
add_subdirectory( src )
add_subdirectory( test )

if( ENABLE_MPAS_DOC )
    add_subdirectory( docs )
endif()

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
