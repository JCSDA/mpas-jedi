# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# MPAS/JEDI bundle
#

cmake_minimum_required( VERSION 3.12 )
project( mpas-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

## ECBuild integration
include(GNUInstallDirs)
set( INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING "Library installation directory" FORCE ) #Force ecbuild to use the default platform lib install location
set( ECBUILD_DEFAULT_BUILD_TYPE Release CACHE STRING "ECBuild default build type" FORCE)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
find_package(ecbuild)
include( ecbuild_bundle )
ecbuild_bundle_initialize()

#TODO: When mpas-bundle becomes a public repo, consider changing the default value of BUNDLE_SKIP_ROPP-UFO to "ON"

option(BUNDLE_SKIP_ROPP-UFO "Don't build ROPP-UFO"  "OFF") # Skip ropp-ufo build unless user passes -DBUNDLE_SKIP_ROPP-UFO=ON
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/jcsda-internal/ropp-test.git" )
ecbuild_bundle( PROJECT crtm      GIT "https://github.com/jcsda/CRTMv3.git" )
ecbuild_bundle( PROJECT rttov     GIT "https://github.com/jcsda-internal/rttov.git" )
ecbuild_bundle( PROJECT oops      GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber     GIT "https://github.com/jcsda-internal/saber.git" )
ecbuild_bundle( PROJECT ioda      GIT "https://github.com/jcsda-internal/ioda.git" )
ecbuild_bundle( PROJECT ufo       GIT "https://github.com/jcsda-internal/ufo.git" )

# ioda and ufo test data
#---------------------------------

ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" )
ecbuild_bundle( PROJECT mpas-jedi-data GIT "https://github.com/JCSDA-internal/mpas-jedi-data.git" )


set(MPAS_CORES atmosphere init_atmosphere CACHE STRING "MPAS-Model: cores to build.")
set(MPAS_DOUBLE_PRECISION "ON" CACHE STRING "MPAS-Model: Use double precision 64-bit Floating point.")

ecbuild_bundle( PROJECT MPAS     GIT "https://github.com/jcsda-internal/MPAS-Model.git" )

option(BUNDLE_SKIP_SACA "Don't build saca"  "ON") # Skip saca build unless user passes -DBUNDLE_SKIP_SACA=OFF

ecbuild_bundle( PROJECT mpas-jedi GIT "https://github.com/jcsda-internal/mpas-jedi")

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
