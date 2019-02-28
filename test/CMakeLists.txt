list( APPEND mpas_test_input
  testinput/3dvar.yaml
  testinput/3dvar_amsua.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/3dvar_bumpcov.existingb.yaml
  testinput/3dvar_bumpcov.parametersbump.yaml
  testinput/3dfgat.yaml
  testinput/3denvar_bump.yaml
  testinput/4denvar_bump.yaml
  testinput/4denvar_ID.yaml
  testinput/3dhybrid_bump.yaml
  testinput/3dhybrid_bump_bump.yaml
  testinput/forecast.yaml
  testinput/hofx.yaml
  testinput/enshofx.yaml
  testinput/dirac.bumploc.yaml
  testinput/dirac.bumpcov.yaml
  testinput/geometry.yaml
  testinput/state.yaml
  testinput/increment.yaml
  testinput/model.yaml
  testinput/errorcov.yaml
  testinput/parametersbump.cov.yaml
  testinput/parametersbump.loc.yaml
  testinput/static_binit.yaml
  testinput/variablechange.yaml
  compare.sh
)

list( APPEND mpas_testoutput
  testoutput/3dvar.test
  testoutput/3dvar_amsua.test
  testoutput/3dvar_bumpcov.test
  testoutput/3dvar_bumpcov.parametersbump.test
  testoutput/3dfgat.test
  testoutput/3denvar_bump.test
  testoutput/4denvar_bump.test
  testoutput/4denvar_ID.test
  testoutput/3dhybrid_bump.test
  testoutput/3dhybrid_bump_bump.test
  testoutput/forecast.test
  testoutput/hofx.test
  testoutput/enshofx.test
  testoutput/dirac.bumploc.test
  testoutput/dirac.bumpcov.test
  testoutput/parametersbump.cov.test
  testoutput/parametersbump.loc.test
  testoutput/staticbinit.test
  testoutput/variablechange.test
)

list( APPEND mpas_graphics
  graphics/plot_obs_nc_loc.py
  graphics/plot_diag_omaomb_profile.py
  graphics/plot_diag_omaomb_timeserial.py
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

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

# Test interface classes with MPAS

ecbuild_add_test( TARGET  test_mpas_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcov.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_variablechange
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/variablechange.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_hofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_hofx.x testinput/hofx.yaml"
                       testoutput/hofx.test
                  MPI     1
                  DEPENDS mpas_hofx.x )

ecbuild_add_test( TARGET  test_mpas_enshofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_enshofx.x testinput/enshofx.yaml"
                       testoutput/enshofx.test
                  MPI     1
                  DEPENDS mpas_enshofx.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumploc.yaml"
                       testoutput/dirac.bumploc.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumpcov
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumpcov.yaml"
                       testoutput/dirac.bumpcov.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_parameters_bump_cov
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_parameters.x testinput/parametersbump.cov.yaml"
                       testoutput/parametersbump.cov.test
                  MPI     1
                  DEPENDS mpas_parameters.x )

ecbuild_add_test( TARGET test_mpas_parameters_bump_loc
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_parameters.x testinput/parametersbump.loc.yaml"
                       testoutput/parametersbump.loc.test
                  MPI     1
                  DEPENDS mpas_parameters.x )

ecbuild_add_test( TARGET test_mpas_forecast
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_forecast.x testinput/forecast.yaml"
                       testoutput/forecast.test
                  MPI     1
                  DEPENDS mpas_forecast.x )

ecbuild_add_test( TARGET test_mpas_3dvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar.yaml"
                       testoutput/3dvar.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dvar_amsua
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar_amsua.yaml"
                       testoutput/3dvar_amsua.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dvar_bumpcov
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar_bumpcov.yaml"
                       testoutput/3dvar_bumpcov.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dvar_bumpcov_parametersbump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar_bumpcov.parametersbump.yaml"
                       testoutput/3dvar_bumpcov.parametersbump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3denvar_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3denvar_bump.yaml"
                       testoutput/3denvar_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dhybrid_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dhybrid_bump.yaml"
                       testoutput/3dhybrid_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dhybrid_bump_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dhybrid_bump_bump.yaml"
                       testoutput/3dhybrid_bump_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dfgat
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dfgat.yaml"
                       testoutput/3dfgat.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_4denvar_ID
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/4denvar_ID.yaml"
                       testoutput/4denvar_ID.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_4denvar_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/4denvar_bump.yaml"
                       testoutput/4denvar_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET  test_mpas_staticbinit
                  TYPE    EXE
                  SOURCES ../mains/mpasStaticBInit.cc
                  ARGS    "testinput/static_binit.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET test_mpas_3dvar_bumpcov_existing_b
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar_bumpcov.existingb.yaml"
                       testoutput/3dvar_bumpcov.test
                  MPI     1
                  DEPENDS mpas_variational.x )
