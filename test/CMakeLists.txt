list( APPEND mpas_test_input
  testinput/3dvar.json
  testinput/3dvar.bumpcov.json
  testinput/3dvar.bumpcov.existingb.json
  testinput/3dfgat.json
  testinput/4denvar_bump.json
  testinput/3dhybrid_bump.json
  testinput/forecast.json
  testinput/hofx.json
  testinput/enshofx.json
  testinput/dirac.bumploc.json
  testinput/dirac.bumpcov.json
  testinput/geometry.json
  testinput/state.json
  testinput/increment.json
  testinput/model.json
  testinput/errorcov.json
  testinput/parametersbump.cov.json
  testinput/parametersbump.loc.json
  testinput/static_binit.json
  testinput/variablechange.json
  compare.sh
)

list( APPEND mpas_testoutput
  testoutput/3dvar.test
  testoutput/3dvar.bumpcov.test
  testoutput/3dfgat.test
  testoutput/4denvar_bump.test
  testoutput/3dhybrid_bump.test
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
                  ARGS    "testinput/geometry.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcov.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_variablechange
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/variablechange.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_hofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_hofx.x testinput/hofx.json"
                       testoutput/hofx.test
                  MPI     1
                  DEPENDS mpas_hofx.x )

ecbuild_add_test( TARGET  test_mpas_enshofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_enshofx.x testinput/enshofx.json"
                       testoutput/enshofx.test
                  MPI     1
                  DEPENDS mpas_enshofx.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumploc.json"
                       testoutput/dirac.bumploc.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumpcov
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumpcov.json"
                       testoutput/dirac.bumpcov.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_parameters_bump_cov
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_parameters.x testinput/parametersbump.cov.json"
                       testoutput/parametersbump.cov.test
                  MPI     1
                  DEPENDS mpas_parameters.x )

ecbuild_add_test( TARGET test_mpas_parameters_bump_loc
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_parameters.x testinput/parametersbump.loc.json"
                       testoutput/parametersbump.loc.test
                  MPI     1
                  DEPENDS mpas_parameters.x )

ecbuild_add_test( TARGET test_mpas_forecast
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_forecast.x testinput/forecast.json"
                       testoutput/forecast.test
                  MPI     1
                  DEPENDS mpas_forecast.x )

ecbuild_add_test( TARGET test_mpas_3dvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar.json"
                       testoutput/3dvar.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dvar_bumpcov
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar.bumpcov.json"
                       testoutput/3dvar.bumpcov.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dfgat
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dfgat.json"
                       testoutput/3dfgat.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_4denvar_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/4denvar_bump.json"
                       testoutput/4denvar_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dhybrid_bump
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dhybrid_bump.json"
                       testoutput/3dhybrid_bump.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET  test_mpas_staticbinit
                  TYPE    EXE
                  SOURCES ../mains/mpasStaticBInit.cc
                  ARGS    "testinput/static_binit.json"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET test_mpas_3dvar_bumpcov_existing_b
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dvar.bumpcov.existingb.json"
                       testoutput/3dvar.bumpcov.test
                  MPI     1
                  DEPENDS mpas_variational.x )
