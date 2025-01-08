# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## cmake functions

#link_files(src_dir dest_dir files..) -Link selected files from source dir to build dir
# Args:
#  src_dir: Full path to source directory
#  dest_dir: Full path to target directory
#  files: A list of file names to link relative to this source dir
macro(link_files src_dir dest_dir)
    foreach(_f IN ITEMS ${ARGN})
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${src_dir}/${_f} ${dest_dir}/${_f})
    endforeach()
endmacro()

#copy_all_files(src_dir dest_dir) - Link all files from source dir to build dir
# Args:
#  src_dir: Full path to source directory
#  dest_dir: Full path to target directory
macro(copy_all_files src_dir dest_dir)
    file(GLOB _all RELATIVE ${src_dir} ${src_dir}/*)
    foreach(_f IN LISTS _all)
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${src_dir}/${_f} ${dest_dir}/)
    endforeach()
endmacro()

#link_all_files(src_dir dest_dir) - Link all files from source dir to build dir
# Args:
#  src_dir: Full path to source directory
#  dest_dir: Full path to target directory
macro(link_all_files src_dir dest_dir)
    file(GLOB _all RELATIVE ${src_dir} ${src_dir}/*)
    foreach(_f IN LISTS _all)
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${src_dir}/${_f} ${dest_dir}/${_f})
    endforeach()
endmacro()

#link_dir(src_dir dest_dir) - Link a full directory
# Args:
#  src_dir: Full path to source directory
#  dest_dir: Full path to target directory link name
macro(link_dir src_dir dest_dir)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${src_dir} ${dest_dir})
endmacro()

# macro to create a symlink from src to dst with just filename
function(CREATE_SYMLINK_FILENAME src dst)
    foreach (FILENAME ${ARGN})
        get_filename_component(filename ${FILENAME} NAME )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${filename} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK_FILENAME)

# Use external jedi-cmake or bundle build
if(DEFINED ENV{jedi_cmake_ROOT})
  include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
else()
  include( ${CMAKE_SOURCE_DIR}/jedicmake/cmake/Functions/git_functions.cmake )
endif()

## end of cmake functions

## test yamls
list( APPEND mpas_testinput
  testinput/3denvar_bumploc.yaml
  testinput/3denvar_dual_resolution.yaml
  testinput/3denvar_2stream_bumploc.yaml
  testinput/3denvar_amsua_allsky.yaml
  testinput/3denvar_amsua_bc.yaml
  testinput/3dfgat.yaml
  testinput/3dfgat_pseudo.yaml
  testinput/3dhybrid_bumpcov_bumploc.yaml
  testinput/3dvar.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/3dvar_bumpcov_ropp.yaml
  testinput/3dvar_bumpcov_nbam.yaml
  testinput/3dvar_bumpcov_rttovcpp.yaml
  testinput/4denvar_bumploc.yaml
  testinput/4denvar_ID.yaml
  testinput/4denvar_VarBC.yaml
  testinput/4denvar_VarBC_nonpar.yaml
  testinput/4dhybrid_bumpcov_bumploc.yaml
  testinput/4dfgat.yaml
  testinput/4dfgat_append_obs.yaml
  testinput/convertstate.yaml
  testinput/converttostructuredgrid_latlon.yaml
  testinput/dirac_bumpcov.yaml
  testinput/dirac_bumploc.yaml
  testinput/dirac_noloc.yaml
  testinput/dirac_spectral_1.yaml
  testinput/dirac_spectral_no_wind.yaml
  testinput/eda_3dhybrid.yaml
  testinput/eda_3dhybrid_1.yaml
  testinput/eda_3dhybrid_2.yaml
  testinput/eda_3dhybrid_3.yaml
  testinput/eda_3dhybrid_4.yaml
  testinput/enshofx.yaml
  testinput/enshofx_1.yaml
  testinput/enshofx_2.yaml
  testinput/enshofx_3.yaml
  testinput/enshofx_4.yaml
  testinput/enshofx_5.yaml
  testinput/errorcovariance.yaml
  testinput/forecast.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/geometry.yaml
  testinput/geometry_iterator_2d.yaml
  testinput/geometry_iterator_3d.yaml
  testinput/getvalues.yaml
  testinput/hofx3d.yaml
  testinput/hofx3d_ropp.yaml
  testinput/hofx3d_nbam.yaml
  testinput/hofx3d_rttovcpp.yaml
  testinput/hofx4d.yaml
  testinput/hofx4d_pseudo.yaml
  testinput/increment.yaml
  testinput/letkf_3dloc.yaml
  testinput/lgetkf.yaml
  testinput/lgetkf_height_vloc.yaml
  testinput/linvarcha.yaml
  testinput/model.yaml
  testinput/obslocalization.yaml
  testinput/obslocalization_vertical.yaml
  testinput/obslocalizations.yaml
  testinput/parameters_bumpcov.yaml
  testinput/parameters_bumploc.yaml
  testinput/process_perts_spectral_no_wind.yaml
  testinput/rtpp.yaml
  testinput/ens_mean_variance.yaml
  testinput/state.yaml
  testinput/unsinterp.yaml
  testinput/obsop_name_map.yaml
)
# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_testinput})

## reference outputs
list( APPEND mpas_testoutput
  testoutput/3denvar_bumploc.ref
  testoutput/3denvar_dual_resolution.ref
  testoutput/3denvar_2stream_bumploc.ref
  testoutput/3denvar_amsua_allsky.ref
  testoutput/3denvar_amsua_bc.ref
  testoutput/3dfgat.ref
  testoutput/3dfgat_pseudo.ref
  testoutput/3dhybrid_bumpcov_bumploc.ref
  testoutput/3dvar.ref
  testoutput/3dvar_bumpcov.ref
  testoutput/3dvar_bumpcov_ropp.ref
  testoutput/3dvar_bumpcov_nbam.ref
  testoutput/3dvar_bumpcov_rttovcpp.ref
  testoutput/4denvar_bumploc.ref
  testoutput/4denvar_ID.ref
  testoutput/4denvar_VarBC.ref
  testoutput/4dhybrid_bumpcov_bumploc.ref
  testoutput/4dfgat.ref
  testoutput/4dfgat_append_obs.ref
  testoutput/convertstate.ref
  testoutput/converttostructuredgrid_latlon.ref
  testoutput/dirac_bumpcov.ref
  testoutput/dirac_bumploc.ref
  testoutput/dirac_noloc.ref
  testoutput/dirac_spectral_1.ref
  testoutput/dirac_spectral_no_wind.ref
  testoutput/eda_3dhybrid.ref
#  testoutput/enshofx.ref
  testoutput/forecast.ref
  testoutput/gen_ens_pert_B.ref
  testoutput/hofx3d.ref
  testoutput/hofx3d_ropp.ref
  testoutput/hofx3d_nbam.ref
  testoutput/hofx3d_rttovcpp.ref
  testoutput/hofx4d.ref
  testoutput/hofx4d_pseudo.ref
  testoutput/letkf_3dloc.ref
  testoutput/lgetkf.ref
  testoutput/lgetkf_height_vloc.ref
  testoutput/parameters_bumpcov.ref
  testoutput/parameters_bumploc.ref
  testoutput/process_perts_spectral_no_wind.ref
  testoutput/rtpp.ref
  testoutput/ens_mean_variance.ref
)
# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_testoutput})

# Create directory for graphics and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/graphics)
link_all_files(${CMAKE_SOURCE_DIR}/mpas-jedi/graphics ${CMAKE_CURRENT_BINARY_DIR}/graphics)

ecbuild_add_resources(TARGET mpas_test_scripts SOURCES_PACK ${mpas_testinput})

############
## Test Data
############
set( DATA_DIR_NAME Data )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME})

##MPAS-Data data table files exported by MPAS-Model core_atmosphere
list( APPEND mpas_static_tables
  CAM_ABS_DATA.DBL
  CAM_AEROPT_DATA.DBL
  COMPATIBILITY
  GENPARM.TBL
  LANDUSE.TBL
  OZONE_DAT.TBL
  OZONE_LAT.TBL
  OZONE_PLEV.TBL
  RRTMG_LW_DATA
  RRTMG_LW_DATA.DBL
  RRTMG_SW_DATA
  RRTMG_SW_DATA.DBL
  SOILPARM.TBL
  VEGPARM.TBL
)
link_files(${MPAS_CORE_ATMOSPHERE_DATADIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_static_tables})

# RTTOV coefficient files
list( APPEND rttov_test_input
  rttov7pred54L/rtcoef_noaa_19_amsua.dat
)

# CRTM coefficient files
list( APPEND crtm_test_input
  AerosolCoeff/Little_Endian/AerosolCoeff.bin
  CloudCoeff/Little_Endian/CloudCoeff.bin
  EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
  EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin
  EmisCoeff/IR_Land/SEcategory/Little_Endian/NPOESS.IRland.EmisCoeff.bin
  EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin
  EmisCoeff/VIS_Ice/SEcategory/Little_Endian/NPOESS.VISice.EmisCoeff.bin
  EmisCoeff/VIS_Land/SEcategory/Little_Endian/NPOESS.VISland.EmisCoeff.bin
  EmisCoeff/VIS_Snow/SEcategory/Little_Endian/NPOESS.VISsnow.EmisCoeff.bin
  EmisCoeff/VIS_Water/SEcategory/Little_Endian/NPOESS.VISwater.EmisCoeff.bin
  EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin
  EmisCoeff/IR_Land/SEcategory/Little_Endian/USGS.IRland.EmisCoeff.bin
  EmisCoeff/VIS_Land/SEcategory/Little_Endian/USGS.VISland.EmisCoeff.bin
  SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_n15.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_n15.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_n18.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_n18.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_aqua.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_aqua.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_metop-a.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_metop-a.TauCoeff.bin
  SpcCoeff/Little_Endian/amsua_metop-b.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/amsua_metop-b.TauCoeff.bin
  SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
  SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
)

# If CRTM is found then we know that UFO must be configured with CRTM as well
# and will have taken care of the data (whether it is provided locally or downloaded);
# in particular it will have set UFO_CRTM_TESTFILES_PATH at parent scope.
if (crtm_FOUND)
  if (NOT UFO_CRTM_TESTFILES_PATH)
    message(FATAL_ERROR "Found crtm but UFO_CRTM_TESTFILES_PATH was not set by ufo")
  elseif (NOT EXISTS ${UFO_CRTM_TESTFILES_PATH})
    message(FATAL_ERROR "Found crtm but ${UFO_CRTM_TESTFILES_PATH}, set by ufo, does not exist")
  endif()
  message(STATUS "Use crtm test data provided by ufo in ${UFO_CRTM_TESTFILES_PATH}")
endif()

# UFO takes care of downloading ufo-data in whatever way specified
# and sets UFO_UFO_TESTFILES_PATH accordingly (includes "/testinput_tier_1")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/ufo)
link_dir(${UFO_UFO_TESTFILES_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/ufo/testinput_tier_1)

# MPAS-JEDI test data - from mpas-jedi-data repo if found, from local directory, or from tarball
if (mpas-jedi-data_FOUND)
  message(STATUS "Use test data from mpas-jedi-data repository")
  set (mpas-jedi_data_testinput_tier_1_local "${CMAKE_SOURCE_DIR}/mpas-jedi-data/testinput_tier_1")
# It's unclear if anyone is using the local data functionality, therefore comment it out
# for now and add in if someone complains / remove if completely if not.
#elseif (DEFINED ENV{MPAS_JEDI_TESTFILES})
#  message(STATUS "Use test data from local directory $ENV{MPAS_JEDI_TESTFILES}")
#  # A bit of guesswork here, I don't know if folks using this option stored it in the same directory structure
#  set (mpas-jedi_data_testinput_tier_1_local "$ENV{MPAS_JEDI_TESTFILES}/mpas-jedi-data/testinput_tier_1")
else()
  set (mpas-jedi_data_testinput_tier_1_local "${CMAKE_SOURCE_DIR}/test-data-release/mpas-jedi/${mpas-jedi_data_tag}/testinput_tier_1")
  if (EXISTS ${mpas-jedi_data_testinput_tier_1_local})
    message(STATUS "Skipping download of test data from ${mpas-jedi_data_download_url}/${mpas-jedi_data_tarball_name} since it already exists")
  else()
    message(STATUS "Download test data from ${mpas-jedi_data_download_url}/${mpas-jedi_data_tarball_name}")
    include(FetchContent)
    FetchContent_Declare(
      mpas-jedi-testdata
      DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/test-data-release
      DOWNLOAD_NO_EXTRACT ON
      URL      ${mpas-jedi_data_download_url}/${mpas-jedi_data_tarball_name}
      URL_HASH MD5=${mpas-jedi_data_tarball_hash}
    )
    FetchContent_MakeAvailable(mpas-jedi-testdata)
    file(ARCHIVE_EXTRACT INPUT ${CMAKE_SOURCE_DIR}/test-data-release/${mpas-jedi_data_tarball_name} DESTINATION ${CMAKE_SOURCE_DIR}/test-data-release)
  endif()
endif()

# link .nc files from mpas-jedi-data to build/mpas-jedi/test/Data
#384km tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/384km)
link_dir(${mpas-jedi_data_testinput_tier_1_local}/384km/init
        ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/384km/init)

#480km
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km)
link_dir(${mpas-jedi_data_testinput_tier_1_local}/480km/bg
         ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km/bg)

#appendable obs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/obs)
link_dir(${mpas-jedi_data_testinput_tier_1_local}/obs/mpasobsappend1
         ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/obs/mpasobsappend1)
link_dir(${mpas-jedi_data_testinput_tier_1_local}/obs/mpasobsappend2
         ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/obs/mpasobsappend2)

#480km_2stream
# .nc files and config files are in the same directory
# therefore .nc files are linked one by one
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km_2stream)
list( APPEND mpas_test_480km_2stream_nc
      mpasout.2018-04-15_00.00.00.nc
      static.nc
      x1.2562.grid.nc
      x1.2562.init.2018-04-14_18.00.00.nc
)
link_files(${mpas-jedi_data_testinput_tier_1_local}/480km_2stream
           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km_2stream
           ${mpas_test_480km_2stream_nc} )

# link streams and namelists from mpas-jedi to build/mpas-jedi/test/Data
#384km tests
list( APPEND mpas_test_384km
      namelist.atmosphere_2018041500
      streams.atmosphere
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/testinput/namelists/384km
           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/384km
           ${mpas_test_384km} )

#480km tests
list( APPEND mpas_test_480km
      namelist.atmosphere_2018041421
      namelist.atmosphere_2018041500
      streams.atmosphere
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/testinput/namelists/480km
           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km
           ${mpas_test_480km} )
# multi_pe_standard_480: list of processor counts for which to run non-enkf multi-PE 480km tests
# options: 2, 4, 6, 8, 12, 16
list(APPEND multi_pe_standard_480 2)
list(APPEND multi_pe_other_480 4 8 16)
list(APPEND mpas_test_partition_480km x1.2562.graph.info)
foreach(THIS_NPE ${multi_pe_standard_480} ${multi_pe_other_480})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_partition_480km x1.2562.graph.info.part.${THIS_NPE})
    endif()
endforeach()
# Link the 480km graph partition files directly into the test directory; required for MPAS-Model
link_files(${CMAKE_CURRENT_SOURCE_DIR}/testinput/namelists/480km
           ${CMAKE_CURRENT_BINARY_DIR}
           ${mpas_test_partition_480km} )

#480km_2stream tests
list( APPEND mpas_test_480km_2stream
      namelist.atmosphere_2018041500
      streams.atmosphere
    )
link_files(${CMAKE_CURRENT_SOURCE_DIR}/testinput/namelists/480km_2stream
           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/480km_2stream
           ${mpas_test_480km_2stream} )

# test config files
list( APPEND mpas_test_config_base
  geovars.yaml
  keptvars.yaml
  stream_list.atmosphere.background
  stream_list.atmosphere.analysis
  stream_list.atmosphere.ensemble
  stream_list.atmosphere.control
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/testinput/namelists
           ${CMAKE_CURRENT_BINARY_DIR}
           ${mpas_test_config_base})

# link crtm coefficient files
if( crtm_FOUND )
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/UFOCoeff)
  create_symlink_filename( ${UFO_CRTM_TESTFILES_PATH}
                           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/UFOCoeff
                           ${crtm_test_input} )
endif ( crtm_FOUND )

# link rttov coefficient files
if( ${rttov_FOUND} )
  create_symlink_filename( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                           ${CMAKE_CURRENT_BINARY_DIR}/Data
                           ${rttov_test_input} )
endif( ${rttov_FOUND} )


# Make directories for output test data
#bump outputs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/bump)
#spectral outputs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/spectral)
#ObsSpaces
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os)
#model states
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states)
#enkf ObsSpaces/mean/perturbations
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/enkf)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/enkf/os)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/enkf/model)

# ensemble member directory names standardized by
# oops/utils/string_utils::swap_name_member
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os/ensemble)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states/ensemble)
set(ensemblemembers mem000;mem001;mem002;mem003;mem004;mem005)
foreach(MEM IN LISTS ensemblemembers)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os/ensemble/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states/ensemble/${MEM})
#  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/enkf/os/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/enkf/model/${MEM})
endforeach(MEM)


## Configure environment variables necessary for running test executables

# The Fortran IO units 101-200 are used by MPAS-Model for unformatted I/O to big_endian binary files
# The MPAS::core::atmosphere target exports the appropriate conversion flags for both gnu and intel,
# Unfortunately the GFortran flag -fconvert only applies when compiling the main program which for
# mpasjedi is a CXX main program which cannot accept the -fconvert flag.  This means the only way to
# instrument the formatting for is via the environment variable GFORTRAN_CONVERT_UNIT which must
# be set for all test executables reading in unformatted data files.
#
# The intel compiler is apparently smarter and the propagated flag `-convert big_endian` from the
# MPAS::core::atmosphere target sets the desired effect on each compiled object file individually,
# so there is no need to use the runtime environment flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(TEST_ENVIRONMENT ENVIRONMENT GFORTRAN_CONVERT_UNIT=big_endian:101-200)
endif()


######################################################################
#
# jedi_unit_test
# ==============
#
# Wrapper to test a JEDI c++ CLASS
#   jedi_unit_test( [ CLASS <name> ]
#                   [ YAMLFILE <filename.yaml> ]
#                   [ NAME <name> ]
#                   [ NPE <number-of-processor-elements> ]
#                   [ TESTTIER <tier-of-test> ])
#
# Options
# -------
#
# CLASS:
#   Name of the C++ class to be tested
#
# YAMLFILE:
#   Name of YAML configuration file (without extension)
#
# NAME: optional
#   A unique test name to use, defaults to be equal to YAMLFILE
#
# NPE: optional
#   Number of processor elements to use in this test
#
# TESTTIER: optional
#   Testing tier for this test
#
######################################################################
function(add_mpasjedi_unit_test)
    set(options )
    set(oneValueArgs CLASS YAMLFILE NAME NPE TESTTIER)
    set(multiValueArgs )
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( NOT DEFINED TEST_TESTTIER )
        set( TEST_TESTTIER "1" )
    endif()

    # If the tier of the test is larger than the current
    # MPASJEDI_TEST_TIER, we exit so it is not added.
    if ( TEST_TESTTIER GREATER MPASJEDI_TEST_TIER)
      return()
    endif()

    if( NOT DEFINED TEST_NAME )
        set( TEST_NAME "${TEST_YAMLFILE}" )
    endif()

    set( MPI_SUFFIX "" )
    if( DEFINED TEST_NPE )
        if( TEST_NPE GREATER 1 )
            set( MPI_SUFFIX "_${TEST_NPE}pe" )
        endif()
    else()
        set( TEST_NPE 1 )
    endif()
    set( TARGET_NAME "${PROJECT_NAME}_${TEST_NAME}${MPI_SUFFIX}" )

    #Configure environment variables
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(_env ENVIRONMENT GFORTRAN_CONVERT_UNIT=big_endian:101-200)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        set(_env ENVIRONMENT F_UFMTENDIAN=little\\$<SEMICOLON>big:101-200)
    endif()

    ecbuild_add_test( TARGET  test_${TARGET_NAME}
                      ${TEST_ENVIRONMENT}
                      SOURCES executables/Test${TEST_CLASS}.cc
                      ARGS    "testinput/${TEST_YAMLFILE}.yaml"
                      MPI     ${TEST_NPE}
                      LIBS    ${PROJECT_NAME})
endfunction()

####################################################################################
#
# add_mpasjedi_application_test
# ==============
#
# Wrapper for JEDI application to
#  (1) run a test
#  (2) then optionally...
#    (a) do nothing (default)
#      OR
#    (b) create a reference output file
#      OR
#    (c) compare new output to a previous reference output
#
#   add_mpasjedi_application_test( [ APPLICATION <name> ]
#                              [ NAME <name> ]
#                              [ RECALIBRATE ]
#                              [ COMPARE ]
#                              [ TOL <float> ]
#                              [ DIF <int> ]
#                              [ NPE <number-of-processor-elements> ]
#                              [ TESTTIER <tier-of-test> ]
#                              [ DEPENDS <previous-test1> [<previous-test2> ...])
#
# Options
# -------
#
# NAME: optional
#   A unique test name. Not needed when APPLICATION will be used for only one
#   test (e.g., forecast). The YAML configuration file must match
#   either NAME, or APPLICATION if NAME is not specified.
#
# APPLICATION:
#   Name of the JEDI application to be tested
#
#
# RECALIBRATE: optional
#   Include RECALIBRATE in order to generate new reference output
#   for future comparisons.
#
# COMPARE: optional
#   Add comparison test in addition to the initial application test.
#   COMPARE enables the comparison of the application output to a
#   reference output.
#
# NPE: optional
#   Number of processor elements to use in this test.
#   Additionally will be used to generate the TARGET_NAME of this test.
#
# N_PARALLEL_STATES: optional
#   Number of parallel model states used in the application.
#     e.g., number of ensemble members in an ensemble application (EnsHofX, EnsForecast, EDA)
#     e.g., number of parallel time windows in 4denvar
#   Additionally will be used to generate the TARGET_NAME of this test.
#
# TESTTIER: optional
#   Testing tier for this test
#
# DEPENDS: optional
#   List of previous application test NAME's on which this test depends.
#   These NAME's and NPE will be used to generate the list of prior
#   TARGET_NAME's (TEST_DEPENDS).
#
####################################################################################
function(add_mpasjedi_application_test)
    set(options RECALIBRATE)
    set(oneValueArgs NAME APPLICATION NPE N_PARALLEL_STATES TESTTIER)
    set(multiValueArgs DEPENDS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( NOT DEFINED TEST_TESTTIER )
        set( TEST_TESTTIER "1" )
    endif()

    # If the tier of the test is larger than the current
    # MPASJEDI_TEST_TIER, we exit so it is not added.
    if ( TEST_TESTTIER GREATER MPASJEDI_TEST_TIER)
      return()
    endif()

    if( NOT DEFINED TEST_NAME )
        set( TEST_NAME "${TEST_APPLICATION}" )
    endif()

    set( MPI_SUFFIX "" )
    if( NOT DEFINED TEST_NPE )
        set( TEST_NPE 1 )
    endif()
    if( NOT DEFINED TEST_N_PARALLEL_STATES )
        set( TEST_N_PARALLEL_STATES 1 )
    endif()
    math(EXPR REMAIN "${TEST_NPE} % ${TEST_N_PARALLEL_STATES}")
    if(NOT ${REMAIN} EQUAL 0)
        message(FATAL_ERROR "NPE % N_PARALLEL_STATES must be zero")
    endif()
    math(EXPR NPE_PER_MEMBER "${TEST_NPE} / ${TEST_N_PARALLEL_STATES}")
    if( ${NPE_PER_MEMBER} GREATER 1 )
        set( MPI_SUFFIX "_${NPE_PER_MEMBER}pe" )
    endif()

    set(TARGET_NAME "${PROJECT_NAME}_${TEST_NAME}${MPI_SUFFIX}")
    set(TARGET_PREFIX "test")
    set(TEST_OUT_FILE "testoutput/${TEST_NAME}.run.ref")
    set(TEST_REFERENCE_FILE "testoutput/${TEST_NAME}.ref")
    if( TEST_RECALIBRATE )
        set( TARGET_PREFIX "setup_ref" )
    endif()
    if( DEFINED TEST_DEPENDS )
        set(templist "")
	foreach(previous_test_NAME ${TEST_DEPENDS})
            list(APPEND templist "${TARGET_PREFIX}_${PROJECT_NAME}_${previous_test_NAME}${MPI_SUFFIX}")
	endforeach(previous_test_NAME)
        set(TEST_DEPENDS "${templist}")
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
			  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_${TEST_APPLICATION}.x
			  ARGS    "testinput/${TEST_NAME}.yaml"
                          MPI     ${TEST_NPE}
                          TEST_DEPENDS "${TEST_DEPENDS}" )
    else()
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
			  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_${TEST_APPLICATION}.x
			  ARGS    "testinput/${TEST_NAME}.yaml"
                          MPI     ${TEST_NPE} )
    endif()
    if( TEST_RECALIBRATE )
        ecbuild_add_test( TARGET  extract_ref_${TARGET_NAME}
                          TYPE    SCRIPT
                          COMMAND cp
                          ARGS    ${TEST_OUT_FILE} ${TEST_REFERENCE_FILE}
                          TEST_DEPENDS ${TARGET_PREFIX}_${TARGET_NAME} )
    endif()
endfunction()

# Default MPASJEDI_TEST_TIER
set(MPASJEDI_TEST_TIER 1)

# Override MPASJEDI_TEST_TIER using environment variable
if(DEFINED ENV{MPASJEDI_TEST_TIER})
    set(MPASJEDI_TEST_TIER "$ENV{MPASJEDI_TEST_TIER}")
endif()

## whether to recalibrate reference files to current build
set(RECALIBRATE_CTEST_REFS "OFF")

if( ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
  set(RECALIBRATE "RECALIBRATE")
else()
  set(RECALIBRATE "")
endif()

#############################
# Individual test definitions
#############################

# UNIT tests of MPAS-JEDI classes
#--------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    # Unit tests for interface classes to PROJECT_NAME
    add_mpasjedi_unit_test( CLASS Geometry        YAMLFILE geometry )
    add_mpasjedi_unit_test( CLASS State           YAMLFILE state )
    add_mpasjedi_unit_test( CLASS Model           YAMLFILE model )
    add_mpasjedi_unit_test( CLASS Increment       YAMLFILE increment )
    add_mpasjedi_unit_test( CLASS ErrorCovariance YAMLFILE errorcovariance )
    add_mpasjedi_unit_test( CLASS LinVarCha       YAMLFILE linvarcha )
    add_mpasjedi_unit_test( CLASS UnstructuredInterpolator YAMLFILE unsinterp NPE 4 )
    add_mpasjedi_unit_test( CLASS GeometryIterator YAMLFILE geometry_iterator_2d NPE 2 )
    add_mpasjedi_unit_test( CLASS GeometryIterator YAMLFILE geometry_iterator_3d NPE 2 )
    add_mpasjedi_unit_test( CLASS GetValues       YAMLFILE getvalues )
    add_mpasjedi_unit_test( CLASS ObsLocalization YAMLFILE obslocalization )
    add_mpasjedi_unit_test( CLASS ObsLocalization YAMLFILE obslocalization_vertical )
    add_mpasjedi_unit_test( CLASS ObsLocalizations YAMLFILE obslocalizations )
endif()

# APPLICATION tests with creation of or comparison to reference output
#---------------------------------------------------------------------
#forecast
add_mpasjedi_application_test(
    APPLICATION forecast
    ${RECALIBRATE})

#hofx3d/hofx4d/enshofx
add_mpasjedi_application_test(
    APPLICATION hofx3d
    ${RECALIBRATE})

add_mpasjedi_application_test(
        NAME hofx3d_nbam
        APPLICATION hofx3d
    ${RECALIBRATE})

if( ${ropp-ufo_FOUND} )
    add_mpasjedi_application_test(
        NAME hofx3d_ropp
        APPLICATION hofx3d
        ${RECALIBRATE})
endif( ${ropp-ufo_FOUND} )

if( ${rttov_FOUND} )
    add_mpasjedi_application_test(
        NAME hofx3d_rttovcpp
        APPLICATION hofx3d
        ${RECALIBRATE})
endif( ${rttov_FOUND} )

add_mpasjedi_application_test(
    NAME hofx4d
    APPLICATION hofx
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME hofx4d_pseudo
    APPLICATION hofx
    ${RECALIBRATE})

### TODO: Get the mpas_enshofx test passing again. Until that time, we remove it.
##add_mpasjedi_application_test(
##    APPLICATION enshofx
##    ${RECALIBRATE})

#convertstate
add_mpasjedi_application_test(
    NAME convertstate
    APPLICATION convertstate
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME converttostructuredgrid_latlon
    APPLICATION converttostructuredgrid
    ${RECALIBRATE})

#parameters
add_mpasjedi_application_test(
    NAME parameters_bumpcov
    APPLICATION error_covariance_toolbox
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME parameters_bumploc
    APPLICATION error_covariance_toolbox
    ${RECALIBRATE})

#gen_ens_pert_B
add_mpasjedi_application_test(
    APPLICATION gen_ens_pert_B
    DEPENDS parameters_bumpcov
    ${RECALIBRATE})

#dirac
add_mpasjedi_application_test(
    NAME dirac_bumpcov
    APPLICATION error_covariance_toolbox
    DEPENDS parameters_bumpcov
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME dirac_bumploc
    APPLICATION error_covariance_toolbox
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME dirac_noloc
    APPLICATION error_covariance_toolbox
    ${RECALIBRATE})

#spectral
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    add_mpasjedi_application_test(
        NAME dirac_spectral_1
        APPLICATION error_covariance_toolbox
        ${RECALIBRATE})
    
    add_mpasjedi_application_test(
        NAME process_perts_spectral_no_wind
        APPLICATION process_perts
        ${RECALIBRATE})

    add_mpasjedi_application_test(
        NAME dirac_spectral_no_wind
        APPLICATION error_covariance_toolbox
        DEPENDS process_perts_spectral_no_wind
        ${RECALIBRATE})
endif()

#variational - 3dvar
add_mpasjedi_application_test(
    NAME 3dvar
    APPLICATION variational
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3dvar_bumpcov
    APPLICATION variational
    DEPENDS parameters_bumpcov
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3dvar_bumpcov_nbam
    APPLICATION variational
    DEPENDS parameters_bumpcov
    ${RECALIBRATE})

if( ${ropp-ufo_FOUND} )
    add_mpasjedi_application_test(
        NAME 3dvar_bumpcov_ropp
        APPLICATION variational
        DEPENDS parameters_bumpcov
        ${RECALIBRATE})
endif( ${ropp-ufo_FOUND} )

if( ${rttov_FOUND} )
    add_mpasjedi_application_test(
        NAME 3dvar_bumpcov_rttovcpp
        APPLICATION variational
        DEPENDS parameters_bumpcov
        ${RECALIBRATE})
endif( ${rttov_FOUND} )

#variational - 3denvar
add_mpasjedi_application_test(
    NAME 3denvar_bumploc
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3denvar_dual_resolution
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3denvar_2stream_bumploc
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3denvar_amsua_allsky
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3denvar_amsua_bc
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

#variational - 3dhybrid
add_mpasjedi_application_test(
    NAME 3dhybrid_bumpcov_bumploc
    APPLICATION variational
    DEPENDS parameters_bumpcov parameters_bumploc
    ${RECALIBRATE})

#variational - 3dfgat
add_mpasjedi_application_test(
    NAME 3dfgat
    APPLICATION variational
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 3dfgat_pseudo
    APPLICATION variational
    ${RECALIBRATE})

#variational - 4denvar
add_mpasjedi_application_test(
    NAME 4denvar_ID
    APPLICATION variational
    NPE 3
    N_PARALLEL_STATES 3
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 4denvar_VarBC
    APPLICATION variational
    NPE 3
    N_PARALLEL_STATES 3
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 4denvar_VarBC_nonpar
    APPLICATION variational
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 4denvar_bumploc
    APPLICATION variational
    NPE 3
    N_PARALLEL_STATES 3
    DEPENDS parameters_bumploc
    ${RECALIBRATE})

#variational - hybrid-4denvar
add_mpasjedi_application_test(
    NAME 4dhybrid_bumpcov_bumploc
    APPLICATION variational
    NPE 3
    N_PARALLEL_STATES 3
    DEPENDS parameters_bumpcov parameters_bumploc
    ${RECALIBRATE})

#variational - 4dfgat
add_mpasjedi_application_test(
    NAME 4dfgat
    APPLICATION variational
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME 4dfgat_append_obs
    APPLICATION variational
    ${RECALIBRATE})

# eda
add_mpasjedi_application_test(
    NAME eda_3dhybrid
    APPLICATION eda
    NPE 4
    N_PARALLEL_STATES 4
    DEPENDS parameters_bumpcov parameters_bumploc
    ${RECALIBRATE})

# rtpp
add_mpasjedi_application_test(
    APPLICATION rtpp
    DEPENDS eda_3dhybrid
    ${RECALIBRATE})

# ens_mean_variance
add_mpasjedi_application_test(
    APPLICATION ens_mean_variance
    DEPENDS eda_3dhybrid
    ${RECALIBRATE})

# enkf (4 and 8 pe)
add_mpasjedi_application_test(
    NAME letkf_3dloc
    NPE 4
    APPLICATION enkf
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME lgetkf
    NPE 4
    APPLICATION enkf
    ${RECALIBRATE})

add_mpasjedi_application_test(
    NAME lgetkf_height_vloc
    APPLICATION enkf
    ${RECALIBRATE})

# subset of APPLICATION tests with comparison to reference output on multiple PE's
#---------------------------------------------------------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    foreach(THIS_NPE ${multi_pe_standard_480})
        if( THIS_NPE GREATER 1 )
	    add_mpasjedi_application_test(
                APPLICATION forecast
                NPE ${THIS_NPE})

## hofx is expensive, hold off list for now
#            add_mpasjedi_application_test(
#                APPLICATION hofx
#                NPE ${THIS_NPE})

            add_mpasjedi_application_test(
                APPLICATION error_covariance_toolbox
		NAME parameters_bumpcov
                NPE ${THIS_NPE})

            add_mpasjedi_application_test(
                APPLICATION error_covariance_toolbox
		NAME parameters_bumploc
                NPE ${THIS_NPE})

            add_mpasjedi_application_test(
                APPLICATION variational
		NAME 3dvar
                NPE ${THIS_NPE})

            add_mpasjedi_application_test(
                APPLICATION variational
		NAME 3dhybrid_bumpcov_bumploc
                NPE ${THIS_NPE}
                DEPENDS parameters_bumpcov parameters_bumploc)
        endif()
    endforeach(THIS_NPE)
endif()
