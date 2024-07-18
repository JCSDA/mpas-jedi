# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# MPAS/JEDI bundle
#

cmake_minimum_required( VERSION 3.14 )
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

# Use external jedi-cmake for findESMF.cmake
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
list( APPEND CMAKE_MODULE_PATH $ENV{jedi_cmake_ROOT}/share/jedicmake/Modules )

# Find external ESMF for mpas-model (optional in theory, but if CI containers have ESMF it becomes mandatory)
find_package(ESMF 8.3.0 MODULE)

#TODO: When mpas-bundle becomes a public repo, consider changing the default value of BUNDLE_SKIP_ROPP-UFO to "ON"

option(BUNDLE_SKIP_ROPP-UFO "Don't build ROPP-UFO"  "OFF") # Skip ropp-ufo build unless user passes -DBUNDLE_SKIP_ROPP-UFO=ON
# Enabling the ioda-data and ufo-data repos but then not cloning them
# when they are not needed saves us from unnecessary downloads of the data tarballs
option(ENABLE_IODA_DATA "Obtain ufo test data from ioda-data repository (vs tarball)" ON)
option(ENABLE_UFO_DATA "Obtain ufo test data from ufo-data repository (vs tarball)" ON)

ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/jcsda-internal/ropp-test.git" )
ecbuild_bundle( PROJECT crtm      GIT "https://github.com/jcsda/CRTMv3.git" )
ecbuild_bundle( PROJECT rttov     GIT "https://github.com/jcsda-internal/rttov.git" )
ecbuild_bundle( PROJECT oops      GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber     GIT "https://github.com/jcsda-internal/saber.git" )
#ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/jcsda-internal/ioda-data.git" )
ecbuild_bundle( PROJECT ioda      GIT "https://github.com/jcsda-internal/ioda.git" )
#ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/jcsda-internal/ufo-data.git" )
ecbuild_bundle( PROJECT ufo       GIT "https://github.com/jcsda-internal/ufo.git" )

# mpas-jedi test data
#---------------------------------
option(ENABLE_MPAS_JEDI_DATA "Obtain mpas-jedi test data from mpas-data repository (vs tarball)" ON)
ecbuild_bundle( PROJECT mpas-jedi-data GIT "https://github.com/JCSDA-internal/mpas-jedi-data.git" )

set(MPAS_CORES atmosphere init_atmosphere CACHE STRING "MPAS-Model: cores to build.")
set(MPAS_DOUBLE_PRECISION "ON" CACHE STRING "MPAS-Model: Use double precision 64-bit Floating point.")

ecbuild_bundle( PROJECT MPAS     GIT "https://github.com/jcsda-internal/MPAS-Model.git" )
ecbuild_bundle( PROJECT mpas-jedi GIT "https://github.com/jcsda-internal/mpas-jedi")

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
