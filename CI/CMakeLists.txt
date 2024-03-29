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

ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" )

# ecmwf repos are already available in CI containers. Disabled here
#option(BUNDLE_SKIP_ECKIT "Don't build eckit" "OFF" ) # Build eckit unless user passes -DBUNDLE_SKIP_ECKIT=ON
#ecbuild_bundle( PROJECT eckit     GIT "https://github.com/ecmwf/eckit.git" )
#ecbuild_bundle( PROJECT fckit     GIT "https://github.com/ecmwf/fckit.git" )
#ecbuild_bundle( PROJECT atlas     GIT "https://github.com/ecmwf/atlas.git" )

#TODO: When mpas-bundle becomes a public repo, consider changing the default value of BUNDLE_SKIP_ROPP-UFO to "ON"

option(BUNDLE_SKIP_ROPP-UFO "Don't build ROPP-UFO"  "OFF") # Skip ropp-ufo build unless user passes -DBUNDLE_SKIP_ROPP-UFO=ON
ecbuild_bundle( PROJECT ropp-ufo  GIT "https://github.com/jcsda-internal/ropp-test.git" )
option(BUNDLE_SKIP_GSW "Don't build GSW"       "ON") # Skip gsw build unless user passes -DBUNDLE_SKIP_GSW=OFF
ecbuild_bundle( PROJECT gsw       GIT "https://github.com/jcsda-internal/GSW-Fortran.git" )
ecbuild_bundle( PROJECT crtm      GIT "https://github.com/jcsda-internal/crtm.git" )
ecbuild_bundle( PROJECT rttov     GIT "https://github.com/jcsda-internal/rttov.git" )
ecbuild_bundle( PROJECT oops      GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber     GIT "https://github.com/jcsda-internal/saber.git" )
ecbuild_bundle( PROJECT ioda      GIT "https://github.com/jcsda-internal/ioda.git" )
ecbuild_bundle( PROJECT ufo       GIT "https://github.com/jcsda-internal/ufo.git" )

# ioda and ufo test data
#---------------------------------

# If IODA branch is being built set GIT_BRANCH_FUNC to IODA's current branch.
# If a tagged version of IODA is being built set GIT_TAG_FUNC to ioda's current tag. In this case,
# IODA test files will be download from UCAR DASH and ioda-data repo will not be cloned.
# When LOCAL_PATH_JEDI_TESTFILES is set to the directory of IODA test files stored
# in a local directory, ioda-data repo will not be cloned

#find_branch_name(REPO_DIR_NAME ioda)
# When LOCAL_PATH_JEDI_TESTFILES is set to the directory of IODA test files stored
# in a local directory, ioda-data repo will not be cloned
#if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
#  ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" BRANCH develop UPDATE )
#endif()

# same procedure for ufo-data
find_branch_name(REPO_DIR_NAME ufo)
if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
  ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" BRANCH develop UPDATE )
endif()

# same procedure for mpas-jedi-data
find_branch_name(REPO_DIR_NAME mpas-jedi)
if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
  ecbuild_bundle( PROJECT mpas-jedi-data GIT "https://github.com/JCSDA-internal/mpas-jedi-data.git" )
endif()


set(MPAS_CORES atmosphere init_atmosphere CACHE STRING "MPAS-Model: cores to build.")
set(MPAS_DOUBLE_PRECISION "ON" CACHE STRING "MPAS-Model: Use double precision 64-bit Floating point.")

ecbuild_bundle( PROJECT MPAS     GIT "https://github.com/jcsda-internal/MPAS-Model.git" )
ecbuild_bundle( PROJECT mpas-jedi GIT "https://github.com/jcsda-internal/mpas-jedi")

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
