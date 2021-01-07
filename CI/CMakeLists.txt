# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# MPAS/JEDI bundle
#

cmake_minimum_required( VERSION 3.12 )
project( mpas-bundle LANGUAGES C CXX Fortran )

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

ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" )

option(BUNDLE_SKIP_ECKIT "Don't build eckit" "ON" ) # Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
ecbuild_bundle( PROJECT eckit     GIT "https://github.com/jcsda-internal/eckit.git" )
ecbuild_bundle( PROJECT fckit     GIT "https://github.com/jcsda-internal/fckit.git" )
ecbuild_bundle( PROJECT atlas     GIT "https://github.com/jcsda-internal/atlas.git" )

#TODO: When mpas-bundle becomes a public repo, consider changing the default value of BUNDLE_SKIP_ROPP-UFO to "ON"
option(BUNDLE_SKIP_ROPP-UFO "Don't build ROPP-UFO"  "OFF") # Skip ropp-ufo build unless user passes -DBUNDLE_SKIP_ROPP-UFO=OFF
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/jcsda-internal/ropp-test.git" )
option(BUNDLE_SKIP_GSW "Don't build GSW"       "ON") # Skip gsw build unless user passes -DBUNDLE_SKIP_GSW=OFF
ecbuild_bundle( PROJECT gsw       GIT "https://github.com/jcsda-internal/GSW-Fortran.git" )
ecbuild_bundle( PROJECT crtm      GIT "https://github.com/jcsda-internal/crtm.git" )

ecbuild_bundle( PROJECT oops      GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber     GIT "https://github.com/jcsda-internal/saber.git" )
ecbuild_bundle( PROJECT ioda      GIT "https://github.com/jcsda-internal/ioda.git" )
ecbuild_bundle( PROJECT ufo       GIT "https://github.com/jcsda-internal/ufo.git" )

set(MPAS_CORES atmosphere init_atmosphere CACHE STRING "MPAS-Model: cores to build.")
set(MPAS_DOUBLE_PRECISION "ON" CACHE STRING "MPAS-Model: Use double precision 64-bit Floating point.")

ecbuild_bundle( PROJECT MPAS-Model GIT "https://github.com/jcsda-internal/MPAS-Model.git" )
ecbuild_bundle( PROJECT mpasjedi GIT "https://github.com/jcsda-internal/mpas-jedi")

ecbuild_bundle_finalize()