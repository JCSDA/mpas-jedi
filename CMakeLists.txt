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
include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )
option( ENABLE_MPAS_DOC "Build MPAS-JEDI documentation" OFF )

## Dependencies
find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS CXX Fortran )
find_package( Boost REQUIRED )

#OpenMP dependent packages
if(OpenMP_FOUND)
    find_package( MPAS 7.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere OpenMP)
    find_package( atlas 0.23.0 REQUIRED COMPONENTS OMP OMP_Fortran)
else()
    find_package( MPAS 7.0 REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere)
    find_package( atlas 0.23.0 REQUIRED )
endif()

find_package( rttov 12.1.0 QUIET )
find_package( oops 1.0.0 REQUIRED )
find_package( saber 1.0.0 REQUIRED )
find_package( ioda 1.0.0 REQUIRED )
find_package( ufo 1.0.0 REQUIRED )

## Sources
add_subdirectory( src )
add_subdirectory( test )

if( ENABLE_MPAS_DOC )
    add_subdirectory( docs )
endif()

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
