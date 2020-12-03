# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

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

## test yamls
list( APPEND mpas_testinput
  testinput/3denvar_bumploc.yaml
  testinput/3dfgat.yaml
  testinput/3dhybrid_bumpcov_bumploc.yaml
  testinput/3dvar.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/4denvar_bumploc.yaml
  testinput/4denvar_ID.yaml
  testinput/dirac_bumpcov.yaml
  testinput/dirac_bumploc.yaml
  testinput/dirac_noloc.yaml
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
  testinput/hofx.yaml
  testinput/hofx_nomodel.yaml
  testinput/increment.yaml
  testinput/linvarcha.yaml
  testinput/model.yaml
  testinput/parameters_bumpcov.yaml
  testinput/parameters_bumploc.yaml
  testinput/state.yaml
  testinput/getvalues.yaml
  testinput/lineargetvalues.yaml
)
# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_testinput})

## reference outputs
list( APPEND mpas_testoutput
  extract_ref.sh
  testoutput/3denvar_bumploc.ref
  testoutput/3dfgat.ref
  testoutput/3dhybrid_bumpcov_bumploc.ref
  testoutput/3dvar.ref
  testoutput/3dvar_bumpcov.ref
  testoutput/4denvar_bumploc.ref
  testoutput/4denvar_ID.ref
  testoutput/dirac_bumpcov.ref
  testoutput/dirac_bumploc.ref
  testoutput/dirac_noloc.ref
  testoutput/eda_3dhybrid.ref
#  testoutput/eda_3denvar.ref
#  testoutput/eda_3dvar.ref
#  testoutput/enshofx.ref
  testoutput/forecast.ref
  testoutput/gen_ens_pert_B.ref
  testoutput/hofx.ref
  testoutput/hofx_nomodel.ref
  testoutput/parameters_bumpcov.ref
  testoutput/parameters_bumploc.ref
)
# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_testoutput})

## graphical scripts
# TODO(JJG): re-organize
list( APPEND mpas_graphics
  graphics/plot_obs_nc_loc.py
  graphics/plot_diag_omaomb.py
  graphics/plot_diag_omaomb_timeserial.py
  graphics/plot_cost_grad.py
  graphics/plot_BUMP_diag.py
  graphics/write_diagnostic_stats.py
  graphics/plot_utils.py
  graphics/plot_stats_timeseries.py
  graphics/plot_inc.py
  graphics/basic_plot_functions.py
)

# This line copies all files to binary dir
#file( COPY ${mpas_test_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# Create directory for graphics and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/graphics)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_graphics})

ecbuild_add_resources(TARGET mpas_test_scripts SOURCES_PACK ${mpas_testinput})

############
## Test Data
############
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#MPAS Model generated namelist, streams, stream_list, and data table files exported by MPAS-Model core_atmosphere
copy_all_files(${MPAS_CORE_ATMOSPHERE_DATADIR} ${CMAKE_CURRENT_BINARY_DIR})


# MPAS static stream lists (modified from model-generated streams)
list( APPEND mpas_test_config_base
  stream_list.atmosphere.output
  stream_list.atmosphere.surface
  stream_list.atmosphere.diagnostics
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/Data  ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_config_base})

## Test files specific to the 480km mesh
# Add or remove values from multi_pe_480 to select multi-PE tests
# 480km options: 2, 4, 6, 8, 12, 16
# streams.atmosphere is modified from the model-generated version for lower resolution
list(APPEND multi_pe_480 2 )
# namelists/streams
list( APPEND mpas_test_config_480
  namelist.atmosphere_2018041421
  namelist.atmosphere_2018041500
  streams.atmosphere
)
# static mesh descriptions
list( APPEND mpas_test_config_480
  x1.2562.graph.info
)
foreach(THIS_NPE ${multi_pe_480})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_config_480 x1.2562.graph.info.part.${THIS_NPE})
    endif()
endforeach()
link_files(${CMAKE_CURRENT_SOURCE_DIR}/Data/480km ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_config_480})

# backgrounds
list( APPEND mpas_test_bg_480
  restart.2018-04-14_21.00.00.nc
  restart.2018-04-15_00.00.00.nc
  restart.2018-04-15_03.00.00.nc
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/Data/480km/bg ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_bg_480})

# ensemble
set( BASE_ENSDIR_480 480km/bg/ensemble )
set(ensmembers mem01;mem02;mem03;mem04;mem05)
list( APPEND mpas_test_ensemble_dates
  2018-04-14_21.00.00
  2018-04-15_00.00.00
  2018-04-15_03.00.00
)
foreach(MEM IN LISTS ensmembers)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/${BASE_ENSDIR_480}/${MEM})
  foreach(DATE ${mpas_test_ensemble_dates})
    set( FILENAME Data/${BASE_ENSDIR_480}/${MEM}/x1.2562.init.${DATE}.nc )
    link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${FILENAME})
  endforeach()
endforeach()

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

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/UFOCoeff)

## Create tests to download IODA, UFO, and CRTM test files

# Test data filenames
set(ioda_test_data ioda_testinput_tier_1.tar.gz)
set(crtm_test_data crtm_coefficients.tar.gz)

# Determine the branch name
# Deafult is to use the name corresponding to the currently checked out branch or tag
# User can override through the MPAS_TESTFILES_BRANCH environment variable

###set( MPAS_TESTFILES_BRANCH "$ENV{MPAS_TESTFILES_BRANCH}")
if (DEFINED ENV{MPAS_TESTFILES_BRANCH})
  message(STATUS "Branch name provided by user: ${MPAS_TESTFILES_BRANCH}")
  set(GIT_BRANCH_MPAS "ENV${MPAS_TESTFILES_BRANCH}")
else()
  # Get the current git branch or tag of mpas repo
  execute_process(
  COMMAND bash "-c" "git symbolic-ref -q --short HEAD || git describe --tags --exact-match"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH_MPAS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "MPAS is in branch: ${GIT_BRANCH_MPAS}")
endif()

# macro for downloading test data
# options:
# - tag 1.0.0: download from UCAR ftp, flat directory
#   Note: when we change the tag, this will have to be updated somehow.
# - anything else: download from AWS, either from branch directory
#   if available, or from develop directory
# input: repo (e.g. ioda or ufo)
#        test_data: test data filename (archive)
# output: output_data_path: path to data on local machine after download
# uses ${GIT_BRANCH_MPAS}
function(ADD_DOWNLOAD_TEST repo test_data output_data_path)
  # for 1.0.0 release always download from UCAR ftp, flat directory
  if(GIT_BRANCH_MPAS STREQUAL "1.0.0")
    set(ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/dataset/147_miesch/file)
    set(download_data_name ${test_data})
    set(data_path ${CMAKE_BINARY_DIR}/test_data/${repo}/1.0.0)
    set(checksum "0")
  # for all other branches/tags download from AWS
  else()
    set(ECBUILD_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com)
    # if branch-specific files exist, download them
    ecbuild_check_urls(NAMES ${repo}/${GIT_BRANCH_MPAS}/${test_data}
                       RESULT MPAS_SPECIFIC_TEST_FILES)
    if(${MPAS_SPECIFIC_TEST_FILES} MATCHES 0)
      set(download_data_name ${repo}/${GIT_BRANCH_MPAS}/${test_data})
      set(data_path ${CMAKE_BINARY_DIR}/test_data/${repo}/${GIT_BRANCH_MPAS})
    # if there are no branch-specific files, download develop files by default
    else()
      set(download_data_name ${repo}/develop/${test_data})
      set(data_path ${CMAKE_BINARY_DIR}/test_data/${repo}/develop)
    endif()
    # only compare checksum if using AWS and not ucar ftp
    set(checksum "1")
  endif()
  message(STATUS "${repo} test data will be downloaded from ${ECBUILD_DOWNLOAD_BASE_URL}/${download_data_name}")
  file(MAKE_DIRECTORY ${data_path})

  # add test for downloading data
  set(testname "mpas_get_${repo}_test_data")
  ecbuild_add_test( TARGET    ${testname}
                    TYPE      SCRIPT
                    COMMAND   ${CMAKE_BINARY_DIR}/bin/mpas_data_downloader.py
                    ARGS      ${download_data_name} ${test_data} ${data_path} ${ECBUILD_DOWNLOAD_BASE_URL} ${checksum})

  # pass output data path back to calling code
  set(${output_data_path} ${data_path} PARENT_SCOPE)
endfunction()

# Currently only use 2.3.0 CRTM coefficients
set( CRTM_COEFFS_BRANCH "2.3.0" )

if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES})
    set(LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}")
endif()

# If local path to testfiles is defined by user don't download
if( DEFINED LOCAL_PATH_JEDI_TESTFILES )
  # Check IODA directory
  if(EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/ioda/${GIT_BRANCH_MPAS})
    set(MPAS_IODA_TESTFILES_PATH ${LOCAL_PATH_JEDI_TESTFILES}/ioda/${GIT_BRANCH_MPAS})
  else()
    set(MPAS_IODA_TESTFILES_PATH ${LOCAL_PATH_JEDI_TESTFILES}/ioda/develop)
  endif()
  message(STATUS "use LOCAL_PATH_JEDI_TESTFILES, ioda files path: ${MPAS_IODA_TESTFILES_PATH}")

  # For CRTM use data with tag 2.3.0
  set(MPAS_CRTM_COEFFS_PATH ${LOCAL_PATH_JEDI_TESTFILES}/crtm/${CRTM_COEFFS_BRANCH})
  message(STATUS "use LOCAL_PATH_JEDI_TESTFILES, crtm coeffs path: ${MPAS_CRTM_COEFFS_PATH}")

# otherwise download the data
else()
  message(STATUS "LOCAL_PATH_JEDI_TESTFILES is not defined, download test files")
  # Create download script for get_ioda_test_data test
  set ( FILENAME mpas_data_downloader.py)
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
  list( APPEND bin_mpas_test_download_scripts_deps ${DEST_FILE} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()

  add_custom_target( bin_mpas_test_download_scripts ALL
      COMMAND chmod +x ${bin_mpas_test_download_scripts_deps}
      DEPENDS ${bin_mpas_test_download_scripts_deps})

  # add tests to download ioda data using ADD_DOWNLOAD_TEST  function
  ADD_DOWNLOAD_TEST( "ioda" ${ioda_test_data} MPAS_IODA_TESTFILES_PATH )

  # download CRTM coefficients
  if( crtm_FOUND )
    set(MPAS_CRTM_COEFFS_PATH ${CMAKE_BINARY_DIR}/test_data/crtm/${CRTM_COEFFS_BRANCH})
    file(MAKE_DIRECTORY ${MPAS_CRTM_COEFFS_PATH})
    # set filenames for downloading: flat directory for release (from UCAR ftp)
    if(GIT_BRANCH_MPAS STREQUAL "1.0.0")
      set(ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/dataset/147_miesch/file)
      set(crtm_download_data_name ${crtm_test_data})
      set(checksum_crtm "0")
    # set filenames for downloading: use directories and branch names for everything else (AWS)
    else()
      set(ECBUILD_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com)
      set(crtm_download_data_name crtm/${CRTM_COEFFS_BRANCH}/${crtm_test_data})
      set(checksum_crtm "1")
    endif()
    ecbuild_add_test( TARGET    mpas_get_crtm_test_data
                      TYPE      SCRIPT
                      COMMAND   ${CMAKE_BINARY_DIR}/bin/mpas_data_downloader.py
                      ARGS      ${crtm_download_data_name} ${crtm_test_data}  ${MPAS_CRTM_COEFFS_PATH} ${ECBUILD_DOWNLOAD_BASE_URL} ${checksum_crtm})

  endif( crtm_FOUND )
endif()

link_dir(${MPAS_IODA_TESTFILES_PATH} ${CMAKE_CURRENT_BINARY_DIR}/Data/ioda)


if( crtm_FOUND )
  create_symlink_filename( ${MPAS_CRTM_COEFFS_PATH}
                           ${CMAKE_CURRENT_BINARY_DIR}/Data/UFOCoeff 
                           ${crtm_test_input} )
endif ( crtm_FOUND )

# Make directories for output test data
#bump outputs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
#ObsSpaces
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/os)
#model states
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/states)
# EDA member directory names standardized by
# oops/utils/string_utils::swap_name_member
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/os/ensemble)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/states/ensemble)
set(ensemblemembers mem001;mem002;mem003;mem004;mem005)
foreach(MEM IN LISTS ensemblemembers)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/os/ensemble/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/states/ensemble/${MEM})
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
    set(TEST_ENVIRONMENT ENVIRONMENT GFORTRAN_CONVERT_UNIT=native\\$<SEMICOLON>big_endian:101-200)
endif()


######################################################################
#
# jedi_unit_test
# ==============
#
# Wrapper to test a JEDI c++ CLASS
#   jedi_unit_test( [ CLASS <name> ]
#                   [ YAMLFILE <filename.yaml> ]
#                   [ ID <id> ]
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
# ID: optional
#   A unique test id to use, defaults to be equal to YAMLFILE
#
# NPE: optional
#   Number of processor elements to use in this test
#
# TESTTIER: optional
#   Testing tier for this test
#
######################################################################
function(jedi_unit_test)
    set(options )
    set(oneValueArgs CLASS YAMLFILE ID NPE TESTTIER)
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

    if( NOT DEFINED TEST_ID )
        set( TEST_ID "${TEST_YAMLFILE}" )
    endif()

    set( MPI_SUFFIX "" )
    if( DEFINED TEST_NPE )
        if( TEST_NPE GREATER 1 )
            set( MPI_SUFFIX "_${TEST_NPE}pe" )
        endif()
    else()
        set( TEST_NPE 1 )
    endif()
    set( TARGET_NAME "${PROJECT_NAME_LOWCASE}_${TEST_ID}${MPI_SUFFIX}" )

    #Configure environment variables
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(_env ENVIRONMENT GFORTRAN_CONVERT_UNIT=native\\$<SEMICOLON>big_endian:101-200)
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
# jedi_app_ref_test_compare
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
#   jedi_app_ref_test_compare( [ APP <name> ]
#                              [ ID <id> ]
#                              [ REFERENCE ]
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
# APP:
#   Name of the JEDI application to be tested
#
# ID: optional
#   A unique test id. Not needed when APP will be used for only one
#   test (e.g., forecast). The YAML configuration file must match
#   either the ID or APP if ID is not specified.
#
# REFERENCE: optional
#   Include REFERENCE in order to generate new reference output
#   for future COMPARE tests. Overrides COMPARE.
#
# COMPARE: optional
#   Add comparison test in addition to the initial application test.
#   COMPARE enables the comparison of the application output to a
#   reference output.
#
# TOL: optional
#   The tolerance for floating-point comparison when executing compare.py
#
# DIF: optional
#   The tolerance for integer comparison when executing compare.py
#   If either TOL or DIF is not defined, compare.sh will be used instead.
#
# NPE: optional
#   Number of processor elements to use in this test.
#   Additionally will be used to generate the TARGET_NAME of this test.
#
# NMEM: optional
#   Number of ensemble members to use in an ensemble application (EnsHofX, EnsForecast, EDA)
#   Additionally will be used to generate the TARGET_NAME of this test.
#
# TESTTIER: optional
#   Testing tier for this test
#
# DEPENDS: optional
#   List of previous application test ID's on which this test depends.
#   These ID's and NPE will be used to generate the list of prior
#   TARGET_NAME's (TEST_DEPENDS).
#
####################################################################################
function(jedi_app_ref_test_compare)
    set(options COMPARE REFERENCE)
    set(oneValueArgs APP ID TOL DIF NPE NMEM TESTTIER)
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

    set( TEST_NAME "${TEST_APP}" )
    if( DEFINED TEST_ID )
        set( TEST_NAME "${TEST_ID}" )
    endif()

    set( MPI_SUFFIX "" )
    if( NOT DEFINED TEST_NPE )
        set( TEST_NPE 1 )
    endif()
    if( NOT DEFINED TEST_NMEM )
        set( TEST_NMEM 1 )
    endif()
    math(EXPR REMAIN "${TEST_NPE} % ${TEST_NMEM}")
    if(NOT ${REMAIN} EQUAL 0)
        message(FATAL_ERROR "NPE % NMEM must be zero")
    endif()
    math(EXPR NPE_PER_MEMBER "${TEST_NPE} / ${TEST_NMEM}")
    if( ${NPE_PER_MEMBER} GREATER 1 )
        set( MPI_SUFFIX "_${NPE_PER_MEMBER}pe" )
    endif()

    set( TARGET_NAME "${PROJECT_NAME_LOWCASE}_${TEST_NAME}${MPI_SUFFIX}" )

    set( TARGET_PREFIX "test" )
    set( TEST_TYPE "run" )
    if( TEST_REFERENCE )
        set( TARGET_PREFIX "setup_ref" )
        set( TEST_TYPE "ref" )
        set( TEST_COMPARE False )
    endif()
    if( DEFINED TEST_DEPENDS )
        set(templist "")
        foreach(previous_test_ID ${TEST_DEPENDS})
           list(APPEND templist "${TARGET_PREFIX}_${PROJECT_NAME_LOWCASE}_${previous_test_ID}${MPI_SUFFIX}")
        endforeach(previous_test_ID)
        set(TEST_DEPENDS "${templist}")
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
                          COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME_LOWCASE}_${TEST_APP}.x
                          ARGS    testinput/${TEST_NAME}.yaml
                                  testoutput/${TEST_NAME}${MPI_SUFFIX}.${TEST_TYPE}
                          MPI     ${TEST_NPE}
                          TEST_DEPENDS "${TEST_DEPENDS}" )
    else()
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
                          COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME_LOWCASE}_${TEST_APP}.x
                          ARGS    testinput/${TEST_NAME}.yaml
                                  testoutput/${TEST_NAME}${MPI_SUFFIX}.${TEST_TYPE}
                          MPI     ${TEST_NPE} )
    endif()
    if( TEST_COMPARE )
        if( DEFINED TEST_TOL AND DEFINED TEST_DIF )
            ecbuild_add_test( TARGET       compare_${TARGET_NAME}
                              TYPE         SCRIPT
                              COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                              ARGS         testoutput/${TEST_NAME}${MPI_SUFFIX}.run
                                           testoutput/${TEST_NAME}.ref
                                           ${TEST_TOL} ${TEST_DIF}
                              TEST_DEPENDS ${TARGET_PREFIX}_${TARGET_NAME} )
        else()
            ecbuild_add_test( TARGET       compare_${TARGET_NAME}
                              TYPE         SCRIPT
                              COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.sh
                              ARGS         testoutput/${TEST_NAME}${MPI_SUFFIX}.run
                                           testoutput/${TEST_NAME}.ref
                              TEST_DEPENDS ${TARGET_PREFIX}_${TARGET_NAME} )
        endif()
    endif()
    if( TEST_REFERENCE )
        ecbuild_add_test( TARGET  extract_ref_${TARGET_NAME}
                          TYPE    SCRIPT
                          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/extract_ref.sh
                          ARGS    testoutput/${TEST_NAME}.ref
                          TEST_DEPENDS ${TARGET_PREFIX}_${TARGET_NAME} )
    endif()
endfunction()

# Default MPASJEDI_TEST_TIER
set(MPASJEDI_TEST_TIER 1)

# Override MPASJEDI_TEST_TIER using environment variable
if(DEFINED ENV{MPASJEDI_TEST_TIER})
    set(MPASJEDI_TEST_TIER "$ENV{MPASJEDI_TEST_TIER}")
endif()

## Settings for comparison to reference output
set(RECALIBRATE_CTEST_REFS "OFF")

# tolerances for compare.py:
set (DEFAULT_FLOAT_TOL 0.00000001) #Default max relative difference
set (FORECAST_FLOAT_TOL 0.005) #Max relative difference for forecast-based ctests
set (INT_TOL 0)   #Max difference

if( ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
  # recalibrate reference files to current build
  set(DEFAULT_REF_TYPE "REFERENCE")
  set(FORECAST_REF_TYPE "REFERENCE")
else()
  # test then compare to reference
  set(DEFAULT_REF_TYPE COMPARE TOL ${DEFAULT_FLOAT_TOL} DIF ${INT_TOL})
  set(FORECAST_REF_TYPE COMPARE TOL ${FORECAST_FLOAT_TOL} DIF ${INT_TOL})
endif()

#############################
# Individual test definitions
#############################

# UNIT tests of MPAS-JEDI classes
#--------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    # Unit tests for interface classes to PROJECT_NAME_LOWCASE
    jedi_unit_test( CLASS Geometry        YAMLFILE geometry )
    jedi_unit_test( CLASS State           YAMLFILE state )
    jedi_unit_test( CLASS Model           YAMLFILE model )
    jedi_unit_test( CLASS Increment       YAMLFILE increment )
    jedi_unit_test( CLASS ErrorCovariance YAMLFILE errorcovariance )
    jedi_unit_test( CLASS LinVarCha       YAMLFILE linvarcha )
    jedi_unit_test( CLASS GetValues       YAMLFILE getvalues )
    jedi_unit_test( CLASS LinearGetValues YAMLFILE lineargetvalues )
endif()

# APPLICATION tests with creation of or comparison to reference output
#---------------------------------------------------------------------
jedi_app_ref_test_compare( APP     forecast                                 ${FORECAST_REF_TYPE} )
jedi_app_ref_test_compare( APP     hofx_nomodel                             ${DEFAULT_REF_TYPE} )
jedi_app_ref_test_compare( APP     hofx                                     ${FORECAST_REF_TYPE} )
### TODO: Get the mpas_enshofx test passing again. Until that time, we remove it.
##jedi_app_ref_test_compare( APP     enshofx                                  ${DEFAULT_REF_TYPE})
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${DEFAULT_REF_TYPE} )
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${DEFAULT_REF_TYPE} )
### TODO: Make the gen_ens_pert_B test produce the same output in Debug and Release. Until that time, we remove it.
#jedi_app_ref_test_compare( APP     gen_ens_pert_B                           ${DEFAULT_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${DEFAULT_REF_TYPE} )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumpcov             ${DEFAULT_REF_TYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumploc             ${DEFAULT_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     dirac       ID dirac_noloc               ${DEFAULT_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar_bumpcov             ${DEFAULT_REF_TYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc           ${DEFAULT_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${DEFAULT_REF_TYPE}
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dfgat                    ${FORECAST_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 4denvar_ID                ${DEFAULT_REF_TYPE} NPE 3 )
jedi_app_ref_test_compare( APP     variational ID 4denvar_bumploc           ${DEFAULT_REF_TYPE} NPE 3
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     eda         ID eda_3dhybrid              ${DEFAULT_REF_TYPE}
                           NPE 4 NMEM 4
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )

# subset of APPLICATION tests with comparison to reference output on multiple PE's
#---------------------------------------------------------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    foreach(THIS_NPE ${multi_pe_480})
        if( THIS_NPE GREATER 1 )
            jedi_app_ref_test_compare( APP     forecast                                 ${FORECAST_REF_TYPE} NPE ${THIS_NPE} )
## hofx is expensive, hold off list for now
#            jedi_app_ref_test_compare( APP     hofx                                     ${FORECAST_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${DEFAULT_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${DEFAULT_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${DEFAULT_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${DEFAULT_REF_TYPE} NPE ${THIS_NPE}
                                       DEPENDS parameters_bumpcov
                                               parameters_bumploc )
        endif()
    endforeach(THIS_NPE)
endif()
