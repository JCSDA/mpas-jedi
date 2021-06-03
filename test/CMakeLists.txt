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

include( ${CMAKE_SOURCE_DIR}/jedicmake/cmake/Functions/git_functions.cmake )

## test yamls
list( APPEND mpas_testinput
  testinput/3denvar_bumploc_bumpinterp.yaml
  testinput/3denvar_bumploc_unsinterp.yaml
  testinput/3denvar_dual_resolution.yaml
  testinput/3denvar_2stream_bumploc_unsinterp.yaml
  testinput/3dfgat.yaml
  testinput/3dhybrid_bumpcov_bumploc.yaml
  testinput/3dvar.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/3dvar_bumpcov_rttovcpp.yaml
  testinput/4denvar_bumploc.yaml
  testinput/4denvar_ID.yaml
  testinput/convertstate_bumpinterp.yaml
  testinput/convertstate_unsinterp.yaml
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
  testinput/hofx3d.yaml
  testinput/hofx3d_rttovcpp.yaml
  testinput/increment.yaml
  testinput/linvarcha.yaml
  testinput/model.yaml
  testinput/parameters_bumpcov.yaml
  testinput/parameters_bumploc.yaml
  testinput/rtpp.yaml
  testinput/state.yaml
  testinput/getvalues_bumpinterp.yaml
  testinput/getvalues_unsinterp.yaml
  testinput/lineargetvalues.yaml
)
# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_testinput})

## reference outputs
list( APPEND mpas_testoutput
  extract_ref.sh
  testoutput/3denvar_bumploc_bumpinterp.ref
  testoutput/3denvar_bumploc_unsinterp.ref
  testoutput/3denvar_dual_resolution.ref
  testoutput/3denvar_2stream_bumploc_unsinterp.ref
  testoutput/3dfgat.ref
  testoutput/3dhybrid_bumpcov_bumploc.ref
  testoutput/3dvar.ref
  testoutput/3dvar_bumpcov.ref
  testoutput/3dvar_bumpcov_rttovcpp.ref
  testoutput/4denvar_bumploc.ref
  testoutput/4denvar_ID.ref
  testoutput/convertstate_bumpinterp.ref
  testoutput/convertstate_unsinterp.ref
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
  testoutput/hofx3d.ref
  testoutput/hofx3d_rttovcpp.ref
  testoutput/parameters_bumpcov.ref
  testoutput/parameters_bumploc.ref
  testoutput/rtpp.ref
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
set( DATA_DIR_NAME Data )
set( BASE_480 480km )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_480})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_480}/bg)

#MPAS Model generated namelist, streams, stream_list, and data table files exported by MPAS-Model core_atmosphere
copy_all_files(${MPAS_CORE_ATMOSPHERE_DATADIR} ${CMAKE_CURRENT_BINARY_DIR})

# MPAS static stream lists (modified from model-generated streams)
list( APPEND mpas_test_config_base
  geovars.yaml
  stream_list.atmosphere.output
  stream_list.atmosphere.surface
  stream_list.atmosphere.diagnostics
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}  ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_config_base})

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
foreach(THIS_NPE ${multi_pe_480})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_partition_480 x1.2562.graph.info.part.${THIS_NPE})
    endif()
endforeach()
# Create 480km links in 480km directory for when we want to specify in yaml
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_480} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_480} ${mpas_test_config_480})
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_480} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_480} ${mpas_test_partition_480})
# Also link just the 480km partition files into the working directory.
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_480} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_partition_480})

# 480km backgrounds
list( APPEND mpas_test_bg_480
  restart.2018-04-14_21.00.00.nc
  restart.2018-04-15_00.00.00.nc
  restart.2018-04-15_03.00.00.nc
)
# Create 480km links in 480km directory for when we want to specify in yaml
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_480}/bg ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_480}/bg ${mpas_test_bg_480})

# ensemble
set( BASE_ENSDIR_480 ${BASE_480}/bg/ensemble )
set(ensmembers mem01;mem02;mem03;mem04;mem05)
list( APPEND mpas_test_ensemble_dates
  2018-04-14_21.00.00
  2018-04-15_00.00.00
  2018-04-15_03.00.00
)
foreach(MEM IN LISTS ensmembers)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_ENSDIR_480}/${MEM})
  foreach(DATE ${mpas_test_ensemble_dates})
    set( FILENAME ${DATA_DIR_NAME}/${BASE_ENSDIR_480}/${MEM}/x1.2562.init.${DATE}.nc )
    link_files(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${FILENAME})
  endforeach()
endforeach()

# 2-stream ctest files
list( APPEND two-stream_480
  analysis.2018-04-15_00.00.00.nc
  mpasout.2018-04-15_00.00.00.nc
  namelist.atmosphere_2018041500
  static.nc
  streams.atmosphere
  x1.2562.grid.nc
)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/480km_2stream)
create_symlink_filename( ${CMAKE_CURRENT_SOURCE_DIR}/Data/480km_2stream
                         ${CMAKE_CURRENT_BINARY_DIR}/Data/480km_2stream
                         ${two-stream_480} )

## Test files specific to the 384km mesh
#  (Used for convertstate test)
# Add or remove values from multi_pe_384 to select multi-PE tests
# 384km options: 2, 4, 6, 8, 12, 16
set( BASE_384 384km )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_384})
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_384}/init)
list(APPEND multi_pe_384 2 )
# namelists/streams
list( APPEND mpas_test_config_384
  namelist.atmosphere_2018041500
  streams.atmosphere
)
# static mesh descriptions
foreach(THIS_NPE ${multi_pe_384})
    if( THIS_NPE GREATER 1 )
        list( APPEND mpas_test_partition_384 x1.4002.graph.info.part.${THIS_NPE})
    endif()
endforeach(THIS_NPE)

# Link all the 384km files into the 384km Data directory
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_384} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_384} ${mpas_test_config_384})
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_384} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_384} ${mpas_test_partition_384})
# Also link just the 384km partition files into the working directory. (No filename conflicts with 480km.)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_384} ${CMAKE_CURRENT_BINARY_DIR} ${mpas_test_partition_384})

# 384km init files
list( APPEND mpas_test_init_384
  x1.4002.init.2018-04-15_00.00.00.nc
)
link_files(${CMAKE_CURRENT_SOURCE_DIR}/${DATA_DIR_NAME}/${BASE_384}/init ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/${BASE_384}/init ${mpas_test_init_384})

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

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/UFOCoeff)

## Create tests to download UFO and CRTM test files


# Currently only use 2.3.0 CRTM coefficients
set( CRTM_COEFFS_BRANCH "2.3.0" )

function(ADD_DOWNLOAD_TEST repo test_files_dirname output_data_path)
  # If REPO branch is being built set GIT_BRANCH_FUNC to REPO's current branch.
  # If a tagged version of REPO is being built set GIT_TAG_FUNC to REPO's current tag.
  find_branch_name(REPO_DIR_NAME ${repo})
  if( DEFINED GIT_BRANCH_FUNC )
    set( REPO_GIT_BRANCH ${GIT_BRANCH_FUNC} )
  elseif( DEFINED GIT_TAG_FUNC )
    set( REPO_GIT_BRANCH ${GIT_TAG_FUNC} )
  endif()

  if ( ${repo} STREQUAL "crtm" )
    set( GIT_TAG_FUNC ${CRTM_COEFFS_BRANCH} )
  endif()

  # When env veriable LOCAL_PATH_JEDI_TESTFILES is set, simply link test files
  # to build directory. get_${repo}_test_data checks the existence of test file directory.
  # ${repo} test data must be stored in
  # ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${BRANCH}/testinput_tier_1
  if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} )
    set( LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}" )
    message( STATUS "use LOCAL_PATH_JEDI_TESTFILES: ${LOCAL_PATH_JEDI_TESTFILES}" )
    # If test data specific to testing branch exists locally use it.
    # If not use test data specific to develop branch.
    if( EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH} )
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH}" )
    else()
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/develop" )
    endif()

    message( STATUS "Test data in ${TESTFILE_DIR_REPO} is linked to build directory" )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
    ${TESTFILE_DIR_REPO} )
    set ( REPO_DATA_DOWNLOADER mpas_data_checker.py )

  # Any tagged version of REPO is being built.
  # Build mpas_get_{repo}_test_data test to download test data from DASH.
  elseif( DEFINED GIT_TAG_FUNC)
    message( STATUS "Tagged version of ${repo} is used" )
    # set ARGS for get_${repo}_test_data
    set( ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/dataset/147_miesch/file )
    set( DIRNAME ${GIT_TAG_FUNC} )
    set( checksum "0" )
    set( TESTFILE_DIR_REPO "${CMAKE_SOURCE_DIR}/test-data-release/${repo}/${DIRNAME}" )

    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${ECBUILD_DOWNLOAD_BASE_URL}
          ${TESTFILE_DIR_REPO}
          ${test_files_dirname}
          ${checksum} )
    # Create test-data-release in source directory
    file( MAKE_DIRECTORY ${TESTFILE_DIR_REPO} )
    set ( REPO_DATA_DOWNLOADER mpas_data_downloader.py )
    message( STATUS "Test data will be downloaded from: ${ECBUILD_DOWNLOAD_BASE_URL}" )

  # Any branch of REPO is being built.
  # repo-data repository is already cloned by bundle/CMakeLists.txt.
  # Link test files in repo-data repo to build directory.
  # get_repo_test_data checks the existence of test file directory
  else()
    set( TESTFILE_DIR_REPO ${CMAKE_SOURCE_DIR}/${repo}-data )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${TESTFILE_DIR_REPO} )
    set (REPO_DATA_DOWNLOADER mpas_data_checker.py )
  endif()

  # add test for downloading data
  ecbuild_add_test( TARGET    mpas_get_${repo}_test_data
                    TYPE      SCRIPT
                    COMMAND   ${CMAKE_BINARY_DIR}/bin/${REPO_DATA_DOWNLOADER}
                    ARGS  ${REPO_DATA_DOWNLOADER_ARGS} )

  # pass output data path back to calling code
  set(${output_data_path} ${TESTFILE_DIR_REPO} PARENT_SCOPE)

  # Create download script for mpas_get_*_test_data test
  set ( FILENAME ${REPO_DATA_DOWNLOADER})
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()


  add_custom_target( bin_mpas_${repo}_test_download_scripts ALL
    COMMAND chmod +x ${DEST_FILE}
    DEPENDS ${DEST_FILE} )

endfunction()

#############################################
# Test data filenames
set(ufo_test_data  ufo_testinput_tier_1.tar.gz)
set(crtm_test_data crtm_coefficients.tar.gz)

ADD_DOWNLOAD_TEST( "ufo" ${ufo_test_data} MPAS_UFO_TESTFILES_PATH )
ADD_DOWNLOAD_TEST( "crtm" ${crtm_test_data} MPAS_CRTM_COEFFS_PATH )



link_dir(${MPAS_UFO_TESTFILES_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/ufo)

if( crtm_FOUND )
  create_symlink_filename( ${MPAS_CRTM_COEFFS_PATH}
                           ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/UFOCoeff
                           ${crtm_test_input} )
endif ( crtm_FOUND )


# Make directories for output test data
#bump outputs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/bump)
#ObsSpaces
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os)
#model states
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states)
# EDA member directory names standardized by
# oops/utils/string_utils::swap_name_member
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os/ensemble)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states/ensemble)
set(ensemblemembers mem000;mem001;mem002;mem003;mem004;mem005)
foreach(MEM IN LISTS ensemblemembers)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/os/ensemble/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DATA_DIR_NAME}/states/ensemble/${MEM})
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
    set( TARGET_NAME "${PROJECT_NAME}_${TEST_ID}${MPI_SUFFIX}" )

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
#                              [ SETUP_REFERENCE ]
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
# SETUP_REFERENCE: optional
#   Include SETUP_REFERENCE in order to generate new reference output
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
#   List of previous application test ID's on which this test depends.
#   These ID's and NPE will be used to generate the list of prior
#   TARGET_NAME's (TEST_DEPENDS).
#
####################################################################################
function(jedi_app_ref_test_compare)
    set(options COMPARE SETUP_REFERENCE)
    set(oneValueArgs APP ID TOL DIF NPE N_PARALLEL_STATES TESTTIER)
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
    set(TEST_LOG_FILE "testoutput/${TEST_NAME}.run")
    set(TEST_REFERENCE_FILE "testoutput/${TEST_NAME}.ref")
    if( TEST_SETUP_REFERENCE )
        set( TARGET_PREFIX "setup_ref" )
    elseif( TEST_COMPARE AND DEFINED TEST_TOL AND DEFINED TEST_DIF )
        ecbuild_add_test( TARGET       compare_${TARGET_NAME}
                          TYPE         SCRIPT
                          COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                          ARGS         ${TEST_LOG_FILE}
                                       ${TEST_REFERENCE_FILE}
                                       ${TEST_TOL} ${TEST_DIF}
                          TEST_DEPENDS ${TARGET_PREFIX}_${TARGET_NAME} )
    endif()
    set(TEST_ARGS "testinput/${TEST_NAME}.yaml")
    if( DEFINED TEST_TOL AND DEFINED TEST_DIF )
      list(APPEND TEST_ARGS "${TEST_LOG_FILE}")
    endif()
    if( DEFINED TEST_DEPENDS )
        set(templist "")
        foreach(previous_test_ID ${TEST_DEPENDS})
           list(APPEND templist "${TARGET_PREFIX}_${PROJECT_NAME}_${previous_test_ID}${MPI_SUFFIX}")
        endforeach(previous_test_ID)
        set(TEST_DEPENDS "${templist}")
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
                          COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_${TEST_APP}.x
                          ARGS    ${TEST_ARGS}
                          MPI     ${TEST_NPE}
                          TEST_DEPENDS "${TEST_DEPENDS}" )
    else()
        ecbuild_add_test( TARGET  ${TARGET_PREFIX}_${TARGET_NAME}
                          ${TEST_ENVIRONMENT}
                          TYPE    EXE
                          COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_${TEST_APP}.x
                          ARGS    ${TEST_ARGS}
                          MPI     ${TEST_NPE} )
    endif()
    if( TEST_SETUP_REFERENCE )
        ecbuild_add_test( TARGET  extract_ref_${TARGET_NAME}
                          TYPE    SCRIPT
                          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/extract_ref.sh
                          ARGS    ${TEST_LOG_FILE}
                                  ${TEST_REFERENCE_FILE}
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
set (PYHX3D_FLOAT_TOL 0.00000001)
set (PYRTTOV_FLOAT_TOL 0.001)
set (PYHX4D_FLOAT_TOL 0.005)
set (PY_INT_TOL 0)   #Max difference

if( ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
  # recalibrate reference files to current build
  set(YAML_REF_TYPE "SETUP_REFERENCE")
  set(PYHX3D_REF_TYPE "SETUP_REFERENCE")
  set(PYRTTOV_REF_TYPE "SETUP_REFERENCE")
  set(PYHX4D_REF_TYPE "SETUP_REFERENCE")
else()
  # test then compare to reference
  set(YAML_REF_TYPE COMPARE)
  set(PYHX3D_REF_TYPE COMPARE TOL ${PYHX3D_FLOAT_TOL} DIF ${PY_INT_TOL})
  set(PYRTTOV_REF_TYPE COMPARE TOL ${PYRTTOV_FLOAT_TOL} DIF ${PY_INT_TOL})
  set(PYHX4D_REF_TYPE COMPARE TOL ${PYHX4D_FLOAT_TOL} DIF ${PY_INT_TOL})
endif()

#############################
# Individual test definitions
#############################

# UNIT tests of MPAS-JEDI classes
#--------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    # Unit tests for interface classes to PROJECT_NAME
    jedi_unit_test( CLASS Geometry        YAMLFILE geometry )
    jedi_unit_test( CLASS State           YAMLFILE state )
    jedi_unit_test( CLASS Model           YAMLFILE model )
    jedi_unit_test( CLASS Increment       YAMLFILE increment )
    jedi_unit_test( CLASS ErrorCovariance YAMLFILE errorcovariance )
    jedi_unit_test( CLASS LinVarCha       YAMLFILE linvarcha )
    jedi_unit_test( CLASS GetValues ID getvalues_bumpinterp YAMLFILE getvalues_bumpinterp )
    jedi_unit_test( CLASS GetValues ID getvalues_unsinterp  YAMLFILE getvalues_unsinterp )
    jedi_unit_test( CLASS LinearGetValues YAMLFILE lineargetvalues )
endif()

# APPLICATION tests with creation of or comparison to reference output
#---------------------------------------------------------------------
jedi_app_ref_test_compare( APP     forecast                                 ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     hofx3d                                   ${PYHX3D_REF_TYPE} )
if( ${rttov_FOUND} )
  jedi_app_ref_test_compare( APP   hofx3d       ID hofx3d_rttovcpp          ${PYRTTOV_REF_TYPE} )
endif( ${rttov_FOUND} )
jedi_app_ref_test_compare( APP     hofx                                     ${PYHX4D_REF_TYPE} )
jedi_app_ref_test_compare( APP     convertstate ID convertstate_bumpinterp  ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     convertstate ID convertstate_unsinterp   ${YAML_REF_TYPE} )
### TODO: Get the mpas_enshofx test passing again. Until that time, we remove it.
##jedi_app_ref_test_compare( APP     enshofx                                  ${YAML_REF_TYPE})
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${YAML_REF_TYPE} )
### TODO: Make the gen_ens_pert_B test produce the same output in Debug and Release. Until that time, we remove it.
#jedi_app_ref_test_compare( APP     gen_ens_pert_B                           ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumpcov             ${YAML_REF_TYPE}
                           DEPENDS parameters_bumpcov )
jedi_app_ref_test_compare( APP     dirac       ID dirac_bumploc             ${YAML_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     dirac       ID dirac_noloc               ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 3dvar_bumpcov             ${YAML_REF_TYPE}
                           DEPENDS parameters_bumpcov )
if( ${rttov_FOUND} )
  jedi_app_ref_test_compare( APP   variational ID 3dvar_bumpcov_rttovcpp    ${YAML_REF_TYPE}
                             DEPENDS parameters_bumpcov )
endif( ${rttov_FOUND} )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc_bumpinterp ${YAML_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3denvar_bumploc_unsinterp ${YAML_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${YAML_REF_TYPE}
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3denvar_dual_resolution   ${YAML_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3denvar_2stream_bumploc_unsinterp ${YAML_REF_TYPE}
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     variational ID 3dfgat                    ${YAML_REF_TYPE} )
jedi_app_ref_test_compare( APP     variational ID 4denvar_ID                ${YAML_REF_TYPE}
                           NPE 3 N_PARALLEL_STATES 3 )
jedi_app_ref_test_compare( APP     variational ID 4denvar_bumploc           ${YAML_REF_TYPE}
                           NPE 3 N_PARALLEL_STATES 3
                           DEPENDS parameters_bumploc )
jedi_app_ref_test_compare( APP     eda         ID eda_3dhybrid              ${YAML_REF_TYPE}
                           NPE 4 N_PARALLEL_STATES 4
                           DEPENDS parameters_bumpcov
                                   parameters_bumploc )
jedi_app_ref_test_compare( APP     rtpp        ID rtpp                      ${DEFAULT_REF_TYPE}
                           DEPENDS eda_3dhybrid )

#####################################################################
# Files for RTTOV tests
#####################################################################
if( ${rttov_FOUND} )
list( APPEND rttov_test_input
rttov7pred54L/rtcoef_noaa_19_amsua.dat
)

# Symlink all RTTOV files
CREATE_SYMLINK_FILENAME( ${rttov_SOURCE_DIR}/rtcoef_rttov12
                         ${CMAKE_CURRENT_BINARY_DIR}/Data
                         ${rttov_test_input} )
endif( ${rttov_FOUND} )

# subset of APPLICATION tests with comparison to reference output on multiple PE's
#---------------------------------------------------------------------------------
if( NOT ${RECALIBRATE_CTEST_REFS} STREQUAL "ON" )
    foreach(THIS_NPE ${multi_pe_480})
        if( THIS_NPE GREATER 1 )
            jedi_app_ref_test_compare( APP     forecast                                 ${YAML_REF_TYPE} NPE ${THIS_NPE} )
## hofx is expensive, hold off list for now
#            jedi_app_ref_test_compare( APP     hofx                                     ${PYHX4D_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumpcov        ${YAML_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     parameters  ID parameters_bumploc        ${YAML_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dvar                     ${YAML_REF_TYPE} NPE ${THIS_NPE} )
            jedi_app_ref_test_compare( APP     variational ID 3dhybrid_bumpcov_bumploc  ${YAML_REF_TYPE} NPE ${THIS_NPE}
                                       DEPENDS parameters_bumpcov
                                               parameters_bumploc )
        endif()
    endforeach(THIS_NPE)
endif()
