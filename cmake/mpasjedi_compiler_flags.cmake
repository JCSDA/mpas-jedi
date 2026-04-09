# (C) Copyright 2026 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


# Set compiler flags for basic build types,
# for compilers where this is not provided by ecbuild.
include(build_type_compiler_flags)

# Set JEDI's common compiler flags
include(jedi_common_compiler_flags)

# Set MPAS-JEDI-specific compiler flags
if(NOT MPAS_DOUBLE_PRECISION)
  add_definitions(-DSINGLE_PRECISION)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL Cray)
  set(CMAKE_CXX_LINK_FLAGS        "-Wl,-Map,loadmap -Wl,-z,muldefs -Ktrap=fp $ENV{CRAYLIBS_X86_64}/btswap.o")
  set(CMAKE_CXX_LINK_EXECUTABLE   "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> <LINK_LIBRARIES> -Wl,-Bdynamic")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  set(CMAKE_Fortran_LINK_FLAGS    "-Wl,-Map,loadmap")
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  ecbuild_add_fortran_flags("-ffpe-trap=invalid,zero,overflow,underflow" BUILD DEBUG)
endif()
