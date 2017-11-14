# (C) Copyright 2017 NCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake)
include( ecbuild_bundle )
project( mpas C CXX Fortran)

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( MPAS_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################

# Boost
set( Boost_MINIMUM_VERSION "1.47" )
find_package( Boost ${Boost_MINIMUM_VERSION} REQUIRED
              COMPONENTS unit_test_framework )
include_directories( ${Boost_INCLUDE_DIR} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.1.0 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# MPAS
include_directories( ${MPAS_INCLUDE} )

################################################################################
# Export package info
################################################################################

list( APPEND MPAS_TPLS oops )

foreach( _tpl ${MPAS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND MPAS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND MPAS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND MPAS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/model )

add_subdirectory( model )
add_subdirectory( test )
add_subdirectory( mains )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME mpas )

# print the summary of the configuration
ecbuild_print_summary()


