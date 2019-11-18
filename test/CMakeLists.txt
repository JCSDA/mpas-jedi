list( APPEND mpas_test_input
  extract_ref.sh
  testinput/3denvar_bumploc.yaml
  testinput/3denvar_bumploc_specified.yaml
  testinput/3dfgat.yaml
  testinput/3dhybrid_bumpcov_bumploc.yaml
  testinput/3dvar.yaml
  testinput/3dvar_amsua.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/3dvar_amsua_bumpcov.yaml
  testinput/4denvar_bumploc.yaml
  testinput/4denvar_ID.yaml
  testinput/dirac_bumpcov.yaml
  testinput/dirac_bumploc.yaml
  testinput/dirac_noloc.yaml
  testinput/enshofx.yaml
  testinput/errorcovariance.yaml
  testinput/forecast.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/geometry.yaml
  testinput/hofx.yaml
  testinput/hofx3d.yaml
  testinput/increment.yaml
  testinput/model.yaml
  testinput/parameters_bumpcov.yaml
  testinput/parameters_bumploc.yaml
  testinput/state.yaml
  testinput/variablechange.yaml
)

list( APPEND mpas_testoutput
  testoutput/3denvar_bumploc.ref
  testoutput/3denvar_bumploc_specified.ref
  testoutput/3dfgat.ref
  testoutput/3dhybrid_bumpcov_bumploc.ref
  testoutput/3dvar.ref
  testoutput/3dvar_amsua.ref
  testoutput/3dvar_amsua_bumpcov.ref
  testoutput/3dvar_bumpcov.ref
  testoutput/4denvar_bumploc.ref
  testoutput/4denvar_ID.ref
  testoutput/dirac_bumpcov.ref
  testoutput/dirac_bumploc.ref
  testoutput/dirac_noloc.ref
#  testoutput/enshofx.ref
  testoutput/forecast.ref
  testoutput/gen_ens_pert_B.ref
  testoutput/hofx.ref
  testoutput/hofx3d.ref
  testoutput/parameters_bumpcov.ref
  testoutput/parameters_bumploc.ref
#  testoutput/staticbinit.ref
#  testoutput/variablechange.ref
)

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

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${mpas_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${mpas_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

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

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

list( APPEND mpas_test_config_base
  Data/stream_list.atmosphere.output
  Data/stream_list.atmosphere.surface
  Data/stream_list.atmosphere.diagnostics
)
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_config_base} )

# Add or remove values from multi_pe_480 to select multi-PE tests
# 480km options: 2, 4, 6, 8, 12, 16
list(APPEND multi_pe_480 2 )

list( APPEND mpas_test_config_480
  Data/480km/namelist.atmosphere_2018041421
  Data/480km/namelist.atmosphere_2018041500
  Data/480km/streams.atmosphere
  Data/480km/x1.2562.graph.info
)
foreach(THIS_NPE ${multi_pe_480})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_config_480
          Data/480km/x1.2562.graph.info.part.${THIS_NPE}
        )
    endif()
endforeach(THIS_NPE)


list( APPEND mpas_test_bg_480
  Data/480km/bg/restart.2018-04-14_21.00.00.nc
  Data/480km/bg/restart.2018-04-15_00.00.00.nc
  Data/480km/bg/restart.2018-04-15_03.00.00.nc
)

list( APPEND mpas_test_ensemble_480
  Data/480km/bg/ensemble/mem01/x1.2562.init.2018-04-14_21.00.00.nc
  Data/480km/bg/ensemble/mem01/x1.2562.init.2018-04-15_00.00.00.nc
  Data/480km/bg/ensemble/mem01/x1.2562.init.2018-04-15_03.00.00.nc
  Data/480km/bg/ensemble/mem02/x1.2562.init.2018-04-14_21.00.00.nc
  Data/480km/bg/ensemble/mem02/x1.2562.init.2018-04-15_00.00.00.nc
  Data/480km/bg/ensemble/mem02/x1.2562.init.2018-04-15_03.00.00.nc
  Data/480km/bg/ensemble/mem03/x1.2562.init.2018-04-14_21.00.00.nc
  Data/480km/bg/ensemble/mem03/x1.2562.init.2018-04-15_00.00.00.nc
  Data/480km/bg/ensemble/mem03/x1.2562.init.2018-04-15_03.00.00.nc
  Data/480km/bg/ensemble/mem04/x1.2562.init.2018-04-14_21.00.00.nc
  Data/480km/bg/ensemble/mem04/x1.2562.init.2018-04-15_00.00.00.nc
  Data/480km/bg/ensemble/mem04/x1.2562.init.2018-04-15_03.00.00.nc
  Data/480km/bg/ensemble/mem05/x1.2562.init.2018-04-14_21.00.00.nc
  Data/480km/bg/ensemble/mem05/x1.2562.init.2018-04-15_00.00.00.nc
  Data/480km/bg/ensemble/mem05/x1.2562.init.2018-04-15_03.00.00.nc
)

# Setup links in test directory for 480km ctests
CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_config_480} )

CREATE_SYMLINK_FILENAME( ${CMAKE_CURRENT_SOURCE_DIR}
                         ${CMAKE_CURRENT_BINARY_DIR}
                         ${mpas_test_bg_480} )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km/bg/ensemble/mem01)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km/bg/ensemble/mem02)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km/bg/ensemble/mem03)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km/bg/ensemble/mem04)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km/bg/ensemble/mem05)
foreach(FILENAME ${mpas_test_ensemble_480})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Add ioda obs test data
list (APPEND ioda_obs_test_data
  atmosphere/aircraft_obs_2018041500_m.nc4
  atmosphere/aircraft_obs_2018041500_s.nc4
  atmosphere/airs_aqua_obs_2018041500_m.nc4
  atmosphere/amsua_n19_obs_2018041500_m.nc4
#  atmosphere/amsua_n19_obs_2018041500_s.nc4
#  atmosphere/atms_npp_obs_2018041500_m.nc4
#  atmosphere/aod_obs_2018041500_m.nc4
#  atmosphere/aod_obs_2018041500_s.nc4
  atmosphere/cris-fsr_npp_obs_2018041500_m.nc4
#  atmosphere/gmi_gpm_obs_2018041500_m.nc4
  atmosphere/gnssro_obs_2018041500_s.nc4
  atmosphere/gnssro_obs_2018041500_m.nc4
#  atmosphere/hirs4_metop-a_obs_2018041500_m.nc4
  atmosphere/iasi_metop-a_obs_2018041500_m.nc4
#  atmosphere/mhs_n19_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_m.nc4
  atmosphere/satwind_obs_2018041500_s.nc4
#  atmosphere/seviri_m08_obs_2018041500_m.nc4
  atmosphere/sfc_obs_2018041500_m.nc4
  atmosphere/sfc_obs_2018041500_s.nc4
#  atmosphere/sndrd1_g15_obs_2018041500_m.nc4
#  atmosphere/sndrd2_g15_obs_2018041500_m.nc4
#  atmosphere/sndrd3_g15_obs_2018041500_m.nc4
#  atmosphere/sndrd4_g15_obs_2018041500_m.nc4
  atmosphere/sondes_obs_2018041500_m.nc4
  atmosphere/sondes_obs_2018041500_s.nc4
  )

CREATE_SYMLINK_FILENAME( ${ioda_SOURCE_DIR}/test/testinput
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${ioda_obs_test_data} )

#####################################################################
# Files for CRTM tests
#####################################################################

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

# Symlink all CRTM files
CREATE_SYMLINK_FILENAME( ${crtm_SOURCE_DIR}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${crtm_test_input} )

#####################################################################

# Make some output directories for test data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis)

######################################################################
#
# jedi_unit_test
# ==============
#
# Wrapper to test a JEDI c++ CLASS
#   jedi_unit_test( [ CLASS <name> ]
#                   [ YAMLFILE <filename.yaml> ]
#                   [ ID <id> ]
#                   [ NPE <number-of-processor-elements> ] )
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
######################################################################
function(jedi_unit_test)
    set(options )
    set(oneValueArgs CLASS YAMLFILE ID NPE)
    set(multiValueArgs )
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

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
#                              [ NPE <number-of-processor-elements> ] )
#                              [ DEPENDS <previous-test1> [<previous-test2> ...]
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
# DEPENDS: optional
#   List of previous application test ID's on which this test depends.
#   These ID's and NPE will be used to generate the list of prior
#   TARGET_NAME's (TEST_DEPENDS).
#
####################################################################################
function(jedi_app_ref_test_compare)
    set(options COMPARE REFERENCE)
    set(oneValueArgs APP ID TOL DIF NPE)
    set(multiValueArgs DEPENDS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    set( TEST_NAME "${TEST_APP}" )
    if( DEFINED TEST_ID )
        set( TEST_NAME "${TEST_ID}" )
    endif() 

    set( MPI_SUFFIX "" )
    if( DEFINED TEST_NPE )
        if( TEST_NPE GREATER 1 )
            set( MPI_SUFFIX "_${TEST_NPE}pe" )
        endif()
    else()
        set( TEST_NPE 1 )
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
    jedi_unit_test( CLASS VariableChange  YAMLFILE variablechange )
endif()

# APPLICATION tests with creation of or comparison to reference output
#---------------------------------------------------------------------
jedi_app_ref_test_compare( APP     forecast                                 ${REFTYPE} )
jedi_app_ref_test_compare( APP     hofx3d                                   ${REFTYPE} )
jedi_app_ref_test_compare( APP     hofx                                     ${REFTYPE} )
### TODO: Get the mpas_enshofx test passing again. Until that time, we remove it.
##jedi_app_ref_test_compare( APP     enshofx                                  ${REFTYPE})
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumpcov             ${REFTYPE} )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumploc             ${REFTYPE} )
jedi_app_ref_test_compare( APP     dirac       ID dirac_noloc               ${REFTYPE} )
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${REFTYPE} )
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${REFTYPE} )
jedi_app_ref_test_compare( APP     gen_ens_pert_B                           ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar_amsua               ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar_bumpcov             ${REFTYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     variational ID 3dvar_amsua_bumpcov       ${REFTYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc           ${REFTYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc_specified ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${REFTYPE}
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dfgat                    ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 4denvar_ID                ${REFTYPE} )
jedi_app_ref_test_compare( APP     variational ID 4denvar_bumploc           ${REFTYPE}
                           DEPENDS parameters_bumploc )

# subset of APPLICATION tests with comparison to reference output on multiple PE's
#---------------------------------------------------------------------------------
if( NOT "${REFTYPE}" STREQUAL "REFERENCE" )
    foreach(THIS_NPE ${multi_pe_480})
        if( THIS_NPE GREATER 1 )
## hofx is expensive, hold off list for now
#            jedi_app_ref_test_compare( APP     hofx                                     ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dvar_amsua               ${COMPARESTR} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${COMPARESTR} NPE ${THIS_NPE}
                                       DEPENDS parameters_bumpcov
                                               parameters_bumploc )
        endif()
    endforeach(THIS_NPE)
endif()

