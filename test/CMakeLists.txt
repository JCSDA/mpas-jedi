# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

## test yamls
list( APPEND mpas_test_input
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
foreach(FILENAME ${mpas_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

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
foreach(FILENAME ${mpas_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

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
foreach(FILENAME ${mpas_graphics})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   mpas_test_scripts
                       SOURCES_PACK
                       ${mpas_test_input}
                     )
############
## Test Data
############
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

## MPAS static stream lists
list( APPEND mpas_test_config_base
  Data/stream_list.atmosphere.output
  Data/stream_list.atmosphere.surface
  Data/stream_list.atmosphere.diagnostics
)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_config_base} )


## Test files specific to the 480km mesh
# Add or remove values from multi_pe_480 to select multi-PE tests
# 480km options: 2, 4, 6, 8, 12, 16
list(APPEND multi_pe_480 2 )
# namelists/streams
list( APPEND mpas_test_config_480
  Data/480km/namelist.atmosphere_2018041421
  Data/480km/namelist.atmosphere_2018041500
  Data/480km/streams.atmosphere
)
# static mesh descriptions
list( APPEND mpas_test_config_480
  Data/480km/x1.2562.graph.info
)
foreach(THIS_NPE ${multi_pe_480})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_config_480
          Data/480km/x1.2562.graph.info.part.${THIS_NPE}
        )
    endif()
endforeach(THIS_NPE)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_config_480} )
# backgrounds
list( APPEND mpas_test_bg_480
  Data/480km/bg/restart.2018-04-14_21.00.00.nc
  Data/480km/bg/restart.2018-04-15_00.00.00.nc
  Data/480km/bg/restart.2018-04-15_03.00.00.nc
)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_bg_480} )
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
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
  endforeach(DATE)
endforeach(MEM)


## IODA observation files
list( APPEND mpas_ioda_test_data
     testinput_tier_1.tar.gz
     )
# Set URL for IODA test files
set( MPAS_IODA_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com)

# Use `LOCAL_PATH_TESTFILES_IODA` to pass path to local test files
if (DEFINED ENV{LOCAL_PATH_TESTFILES_IODA})
    set(LOCAL_PATH_TESTFILES_IODA "$ENV{LOCAL_PATH_TESTFILES_IODA}")
endif()

# If it is not defined then download to build bundle
# in `test_data` directory
# If local path to testfiles is defined don't download

if( NOT DEFINED LOCAL_PATH_TESTFILES_IODA )
  # Get the current git branch of mpas-jedi and use it to
  # download specific data for this branch on S3. If specific
  # data does not exist, it will use default (develop) data

  if( NOT DEFINED MPAS_IODA_TESTFILES_BRANCH)
    # MPAS_IODA_TESTFILES_BRANCH can be used when building
    # to set a specific branch name for testing.

    execute_process(
      COMMAND git rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH_MPAS
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    message( STATUS "Branch name provided by user")
    set( GIT_BRANCH_MPAS ${MPAS_IODA_TESTFILES_BRANCH})
  endif()

   # Check whether the URLs exist or not
  set( ECBUILD_DOWNLOAD_BASE_URL ${MPAS_IODA_DOWNLOAD_BASE_URL}/ioda)
   foreach (mpas_ioda_tar_name ${mpas_ioda_test_data})
       list(APPEND url_names_data ${GIT_BRANCH_MPAS}/${mpas_ioda_tar_name})
   endforeach()
   ecbuild_check_urls(NAMES  ${url_names_data}
                      RESULT MPAS_SPECIFIC_TEST_FILES)

  # Set distant directory
  if(${MPAS_SPECIFIC_TEST_FILES} MATCHES 0)
    # Download and extract new test files (distant directory = git branch)
    set(DIRNAME ${GIT_BRANCH_MPAS})
  else()
    # Download and extract develop test files (distant directory = develop)
    set(DIRNAME "develop")
  endif()
  message(STATUS "Test data will be downloaded from: " ${ECBUILD_DOWNLOAD_BASE_URL}/${DIRNAME})

  set(MPAS_IODA_REP_NAME ioda)
  set(MPAS_IODA_TESTFILES_NAME ${mpas_ioda_test_data})
  set(MPAS_IODA_BRANCH_NAME ${DIRNAME})
  set(MPAS_IODA_TESTFILES_PATH ${CMAKE_BINARY_DIR}/test_data/ioda/${DIRNAME})
  message(STATUS "Save data to " ${MPAS_IODA_TESTFILES_PATH})
  file(MAKE_DIRECTORY ${MPAS_IODA_TESTFILES_PATH})

  # Create download script for get_ioda_test_data test
  set ( FILENAME mpas_ioda_data_downloader.py)
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
  list( APPEND bin_mpas_ioda_test_download_scripts_deps ${DEST_FILE} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()

  add_custom_target( bin_mpas_ioda_test_download_scripts ALL
      COMMAND chmod +x ${bin_mpas_ioda_test_download_scripts_deps}
      DEPENDS ${bin_mpas_ioda_test_download_scripts_deps})

  ecbuild_add_test( TARGET    mpas_get_ioda_test_data
                    TYPE      SCRIPT
                    COMMAND    ${CMAKE_BINARY_DIR}/bin/${FILENAME}
                    ARGS      testoutput/download_test_mpas.log)

else()
  set(MPAS_IODA_TESTFILES_PATH ${LOCAL_PATH_TESTFILES_IODA})
  message (STATUS "LOCAL_PATH_TESTFILES_IODA set to : " ${LOCAL_PATH_TESTFILES_IODA})
endif()

# Link test files
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${MPAS_IODA_TESTFILES_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/ioda )


## CRTM coefficient files
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
#  SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
#  SpcCoeff/Little_Endian/hirs4_metop-a.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/hirs4_metop-a.TauCoeff.bin
  SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
  TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
#  SpcCoeff/Little_Endian/mhs_n19.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/mhs_n19.TauCoeff.bin
#  SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
#  TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
#  SpcCoeff/Little_Endian/sndrD1_g15.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/sndrD1_g15.TauCoeff.bin
#  SpcCoeff/Little_Endian/sndrD2_g15.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/sndrD2_g15.TauCoeff.bin
#  SpcCoeff/Little_Endian/sndrD3_g15.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/sndrD3_g15.TauCoeff.bin
#  SpcCoeff/Little_Endian/sndrD4_g15.SpcCoeff.bin
#  TauCoeff/ODPS/Little_Endian/sndrD4_g15.TauCoeff.bin
#  SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
#  TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
)
CREATE_SYMLINK_FILENAME( ${crtm_SOURCE_DIR}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${crtm_test_input} )

#####################################################################

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


    ecbuild_add_test( TARGET  test_${TARGET_NAME}
                      SOURCES executables/Test${TEST_CLASS}.cc
                      ARGS    "testinput/${TEST_YAMLFILE}.yaml"
                      MPI     ${TEST_NPE}
                      LIBS    ${PROJECT_NAME_LOWCASE} )
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
                          TYPE    EXE
                          COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME_LOWCASE}_${TEST_APP}.x
                          ARGS    testinput/${TEST_NAME}.yaml
                                  testoutput/${TEST_NAME}${MPI_SUFFIX}.${TEST_TYPE}
                          MPI     ${TEST_NPE}
                          TEST_DEPENDS "${TEST_DEPENDS}" )
    else()
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
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

## Settings for comparison to reference output
# tolernaces for compare.py:
set (ctol 0.0) #Max relative difference
set (idif 0)   #Max difference

# use compare.py:
set(COMPARESTR COMPARE TOL ${ctol} DIF ${idif})
# use compare.sh:
#set(COMPARESTR COMPARE)

# Default MPASJEDI_TEST_TIER
set(MPASJEDI_TEST_TIER 1)

# Override MPASJEDI_TEST_TIER using environment variable
if(DEFINED ENV{MPASJEDI_TEST_TIER})
    set(MPASJEDI_TEST_TIER "$ENV{MPASJEDI_TEST_TIER}")
endif()

################
# Select REFTYPE
################
## For APPLICATION tests, choose either
## (a) create reference files (REFTYPE == "REFERENCE"):
#set(REFTYPE "REFERENCE")
##
##      OR
##
## (b) test then compare to reference for all tests:
set(REFTYPE "${COMPARESTR}")

#############################
# Individual test definitions
#############################

# UNIT tests of MPAS-JEDI classes
#--------------------------------
if( NOT "${REFTYPE}" STREQUAL "REFERENCE" )
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
jedi_app_ref_test_compare( APP     forecast                                 ${REFTYPE} )
jedi_app_ref_test_compare( APP     hofx_nomodel                             ${REFTYPE} )
jedi_app_ref_test_compare( APP     hofx                                     ${REFTYPE} )
### TODO: Get the mpas_enshofx test passing again. Until that time, we remove it.
##jedi_app_ref_test_compare( APP     enshofx                                  ${REFTYPE})
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${REFTYPE} )
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${REFTYPE} )
### TODO: Make the gen_ens_pert_B test produce the same output in Debug and Release. Until that time, we remove it.
#jedi_app_ref_test_compare( APP     gen_ens_pert_B                           ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${REFTYPE} )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumpcov             ${REFTYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumploc             ${REFTYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     dirac       ID dirac_noloc               ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar_bumpcov             ${REFTYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc           ${REFTYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${REFTYPE}
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dfgat                    ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 4denvar_ID                ${REFTYPE} NPE 3 )
jedi_app_ref_test_compare( APP     variational ID 4denvar_bumploc           ${REFTYPE} NPE 3
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     eda         ID eda_3dhybrid              ${REFTYPE}
                           NPE 4 NMEM 4
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )

# subset of APPLICATION tests with comparison to reference output on multiple PE's
#---------------------------------------------------------------------------------
if( NOT "${REFTYPE}" STREQUAL "REFERENCE" )
    foreach(THIS_NPE ${multi_pe_480})
        if( THIS_NPE GREATER 1 )
## hofx is expensive, hold off list for now
#            jedi_app_ref_test_compare( APP     hofx                                     ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${COMPARESTR} NPE ${THIS_NPE}
                                       DEPENDS parameters_bumpcov
                                               parameters_bumploc )
        endif()
    endforeach(THIS_NPE)
endif()
