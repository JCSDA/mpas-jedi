#!/bin/csh -f
# set echo
#------------------------------------------------------------------------------#
# Building jedi with a the MPAS model interface
#
# Make sure directory "jedi" has been checked out first with the command:
#   git clone https://github.com/UCAR/oops.git
#
# Make sure to be on the nicas branch with the command:
#   git checkout feature/nicas
#
# Make sure directory "mpas" has been checked out first with the command:
#    git clone https://github.com/UCAR/mpas.git
#
# Make sure directories "mpas", "ufo" and "oops" are at the same level in the top dir
#
# Build will be done in ~/jedi
#
# Prerequies:
# 1. make ufo https://github.com/UCAR/ufo
#    git ioda, crtm, oops in ~/jedi/code
#    make mkufo
#    build_ufo.csh
#
#
#------------------------------------------------------------------------------#
# Define environment
#------------------------------------------------------------------------------#
# This is where the model sources were checked out
    # working with singularity, code sprint Nov 2017
    #OOPS ~/jedi/code/oops
    #PWD ~/jedi/code/mpas 
    #BUILD ~/jedi/build
    setenv MODEL "mpas"
    setenv OOPS "oops" 
    setenv SRC_MODEL $PWD
    echo "SRC_MODEL $SRC_MODEL"
    # MPAS Need to be pre-compiled. If MPAS not needed,set MPAS_INCLUDE=0
    setenv SRC_MPAS "/home/vagrant/jedi_ufo/code/MPAS-Release/src"

# This is where OOPS sources were checked out
    setenv SRC_OOPS  $PWD/../$OOPS
    echo "$SRC_OOPS"

# This is where to build the code
    setenv BUILD     $PWD/../../build
    echo "$BUILD"

#------------------------------------------------------------------------------#
# Compilation module and libraries
#------------------------------------------------------------------------------#
# We do not need so much using pre-set singularity environement
# add for cheyenne
 
# Define lapack path
#setenv LAPACK_LIBRARIES "$LAPACK_PATH/lib64/liblapack.a;$LAPACK_PATH/lib64/libblas.a"
setenv LAPACK_LIBRARIES "$LAPACK_LIBRARIES"

# Need eigen3 library
#setenv EIGEN3_INCLUDE_DIR "" #/glade/p/ral/nsap/jcsda/code/eigen/build

# Need boost library
#setenv BOOST_ROOT /glade/p/ral/nsap/jcsda/code/boost_1_64_0

#Need NETCDF library
#setenv NETCDF_LIBRARIES "${NETCDF}/lib/libnetcdf.a;${NETCDF}/lib/libnetcdff.a"
setenv NETCDF_LIBRARIES "${NETCDF}/lib/libnetcdf.a;${NETCDF}/lib/libnetcdff.a"

setenv MPAS_LIBRARIES "${SRC_MPAS}/libframework.a;${SRC_MPAS}/libdycore.a;${SRC_MPAS}/libops.a"
setenv MPAS_INCLUDE "${SRC_MPAS}/framework;${SRC_MPAS}/operators;${SRC_MPAS}/core_atmosphere"

setenv UFO_LIBRARIES="${BUILD}/ufo"

#------------------------------------------------------------------------------#
# Building oops first setting build_oops=1
#------------------------------------------------------------------------------#
set build_oops=0

if ( $build_oops ) then
# Set path to find ecbuild
  set path = (${path} ${SRC_OOPS}/ecbuild/bin $EIGEN3_INCLUDE_DIR)

# Clean-up directory build/jedi
  rm -rf ${BUILD}/${OOPS}; mkdir -p ${BUILD}/${OOPS}; cd ${BUILD}/${OOPS}

# configure on cheyenne
#  ecbuild --build=debug -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON -DLAPACK_PATH=$LAPACK_PATH -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES -NETCDF_LIBRARIES=${NETCDF_LIBRARIES} -DNETCDF_PATH=${NETCDF} ${SRC_OOPS}

  ecbuild --build=debug -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES ${SRC_OOPS}
# Compile
  make VERBOSE=1 -j4

endif


#------------------------------------------------------------------------------#
# Building model (eg: mpas): if oops is pre-built, set  build_oops=0
#------------------------------------------------------------------------------#

# Set path to find ecbuild
#  set path = (${path} ${SRC_OOPS}/ecbuild/bin )# ${SRC_OOPS}/ecbuild/share/ecbuild/cmake)

# Clean-up directory build/jedi
  rm -rf ${BUILD}/${MODEL}; mkdir ${BUILD}/${MODEL}; cd ${BUILD}/${MODEL}

    echo "SRC_OOPS:  $SRC_OOPS"
    echo "SRC_MODEL: $SRC_MODEL"
    echo "BUILD:     $BUILD"
    echo "OOPS COMPILED: ${BUILD}/${OOPS}"
    echo "MPAS LIBS:     ${MPAS_LIBRARIES}"
    echo "MPAS INCLUDE:  ${MPAS_INCLUDE}"
    echo "UFO PATH:  ${BUILD}/ufo"
# Configure without ufo
#  ecbuild --build=debug -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON -DLAPACK_PATH=$LAPACK_PATH -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES -DNETCDF_LIBRARIES=${NETCDF_LIBRARIES} -DNETCDF_PATH=${NETCDF} -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDE=$MPAS_INCLUDE -DOOPS_PATH=${BUILD}/${OOPS} ${SRC_MODEL}

setenv UFO_PATH "~/jedi/build/ufo"
setenv UFO_INCLUDE "${UFO_PATH}/module"

#ecbuild --build=debug -DUFO_PATH=${UFO_PATH} -DUFO_INCLUDE=${UFO_INCLUDE} -DOOPS_PATH=${BUILD}/${OOPS} -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDE=$MPAS_INCLUDE  ${SRC_MODEL}

#ecbuild --build=debug -DUFO_PATH=${UFO_PATH} -DUFO_INCLUDE=${UFO_INCLUDE} -DUFO_LIBRARIES=${UFO_LIBRARIES} -DOOPS_PATH=${BUILD}/${OOPS} -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDE=$MPAS_INCLUDE  ${SRC_MODEL}

ecbuild --build=debug -DUFO_PATH=${UFO_PATH} -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDE=$MPAS_INCLUDE  ${SRC_MODEL}

#cd ${BUILD}/${MODEL}
#make -j4

exit
# Compile
 make -j4

exit 0

#------------------------------------------------------------------------------#
# The following is needed at run time on Cheyenne
#------------------------------------------------------------------------------#
   unsetenv LD_LIBRARY_PATH module purge
   module load gnu mpich/3.2 cmake/3.7.2 netcdf
   setenv BOOST_ROOT /glade/p/ral/nsap/jcsda/code/boost_1_64_0
   setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${BOOST_ROOT}/stage/lib"
