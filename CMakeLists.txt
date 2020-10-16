# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# MPAS 
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( mpasjedi LANGUAGES C CXX Fortran )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )

## Dependencies
find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS CXX Fortran )
find_package( Boost REQUIRED )

#OpenMP dependent packages
if(OpenMP_FOUND)
    find_package( MPAS REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere OpenMP)
    find_package( atlas REQUIRED COMPONENTS OMP OMP_Fortran)
else()
    find_package( MPAS REQUIRED COMPONENTS DOUBLE_PRECISION core_atmosphere)
    find_package( atlas REQUIRED )
endif()

find_package( oops REQUIRED )
find_package( saber REQUIRED )
find_package( ufo REQUIRED )

## Sources
add_subdirectory( src/mpasjedi )
add_subdirectory( src/mains )
add_subdirectory( test )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
