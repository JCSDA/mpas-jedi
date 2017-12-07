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
# Make sure directory "wrf" has been checked out first with the command:
#    git clone https://github.com/JCSDA/wrf.git
#
# Make sure directories "mpas" and "jedi" are at the same level in the top dir
#
# Build will be done in ../build_jedi
#
#------------------------------------------------------------------------------#
# Define environment
#------------------------------------------------------------------------------#
setenv MODEL "mpas"
setenv SRC_MPAS "/home/vagrant/jedi/code/jedi-bundle/MPAS-Release/src"
setenv MPAS_LIBRARIES "${SRC_MPAS}/libframework.a;${SRC_MPAS}/libdycore.a;${SRC_MPAS}/libops.a"
setenv MPAS_INCLUDES "${SRC_MPAS}/driver;${SRC_MPAS}/framework;${SRC_MPAS}/operators;${SRC_MPAS}/core_atmosphere"
echo "MPAS LIBS:     ${MPAS_LIBRARIES}"
echo "MPAS INCLUDE:  ${MPAS_INCLUDES}"

setenv SRC_MODEL "/home/vagrant/jedi/code/mpas-bundle/${MODEL}"
setenv BUILD "/home/vagrant/jedi/build_test/mpas-bundle/${MODEL}"
echo "BUILD $BUILD"
rm -rf $BUILD
mkdir -p $BUILD
cd $BUILD

#ecbuild --build=debug -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON -DLAPACK_PATH=$LAPACK_PATH -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES -DNETCDF_LIBRARIES=${NETCDF_LIBRARIES} -DNETCDF_PATH=${NETCDF} -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDE=$MPAS_INCLUDE -DOOPS_PATH=${BUILD}/${OOPS} ${SRC_MODEL}
#  ecbuild --build=debug -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES -DOOPS_PATH=${BUILD}/jedi/${OOPS} ${SRC_MODEL}

ecbuild /home/vagrant/jedi/code/mpas-bundle

#ecbuild -DMPAS_LIBRARIES=${MPAS_LIBRARIES} -DMPAS_INCLUDES=$MPAS_INCLUDES /home/vagrant/jedi/code/mpas-bundle
# Compile
 make -j4

exit 0

#------------------------------------------------------------------------------#
# The following is needed at run time:
#------------------------------------------------------------------------------#
   unsetenv LD_LIBRARY_PATH module purge
   module load gnu mpich/3.2 cmake/3.7.2 netcdf
   setenv BOOST_ROOT /glade/p/ral/nsap/jcsda/code/boost_1_64_0
   setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${BOOST_ROOT}/stage/lib"
