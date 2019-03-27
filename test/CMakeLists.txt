list( APPEND mpas_test_input
  testinput/3dvar.yaml
  testinput/3dvar_amsua.yaml
  testinput/3dvar_bumpcov.yaml
  testinput/3dfgat.yaml
  testinput/3denvar_bumploc.yaml
  testinput/3denvar_bumploc_specified.yaml
  testinput/4denvar_bumploc.yaml
  testinput/4denvar_ID.yaml
  testinput/3dhybrid_bumpcov_bumploc.yaml
  testinput/forecast.yaml
  testinput/hofx.yaml
  testinput/hofx3d.yaml
  testinput/enshofx.yaml
  testinput/dirac.bumpcov.yaml
  testinput/dirac.bumploc.yaml
  testinput/dirac.noloc.yaml
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
  testoutput/3dfgat.test
  testoutput/3denvar_bumploc.test
  testoutput/3denvar_bumploc_specified.test
  testoutput/4denvar_bumploc.test
  testoutput/4denvar_ID.test
  testoutput/3dhybrid_bumpcov_bumploc.test
  testoutput/forecast.test
  testoutput/hofx.test
  testoutput/hofx3d.test
  testoutput/enshofx.test
  testoutput/dirac.bumpcov.test
  testoutput/dirac.bumploc.test
  testoutput/dirac.noloc.test
  testoutput/parametersbump.cov.test
  testoutput/parametersbump.loc.test
#  testoutput/staticbinit.test
#  testoutput/variablechange.test
)

list( APPEND mpas_graphics
  graphics/plot_obs_nc_loc.py
  graphics/plot_diag_omaomb.py
  graphics/plot_diag_omaomb_timeserial.py
  graphics/plot_cost_grad.py
  graphics/plot_BUMP_diag.py
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
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_state
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_model
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcov.yaml"
                  MPI     1
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_variablechange
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

ecbuild_add_test( TARGET  test_mpas_hofx3d
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_hofx3d.x testinput/hofx3d.yaml"
                       testoutput/hofx3d.test
                  MPI     1
                  DEPENDS mpas_hofx3d.x )

ecbuild_add_test( TARGET  test_mpas_enshofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_enshofx.x testinput/enshofx.yaml"
                       testoutput/enshofx.test
                  MPI     1
                  DEPENDS mpas_enshofx.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumpcov
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumpcov.yaml"
                       testoutput/dirac.bumpcov.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_dirac_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.bumploc.yaml"
                       testoutput/dirac.bumploc.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_dirac_noloc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.noloc.yaml"
                       testoutput/dirac.noloc.test
                  MPI     1
                  DEPENDS mpas_dirac.x )

ecbuild_add_test( TARGET test_mpas_parameters_bumpcov
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_parameters.x testinput/parametersbump.cov.yaml"
                       testoutput/parametersbump.cov.test
                  MPI     1
                  DEPENDS mpas_parameters.x )

ecbuild_add_test( TARGET test_mpas_parameters_bumploc
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

ecbuild_add_test( TARGET test_mpas_3denvar_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3denvar_bumploc.yaml"
                       testoutput/3denvar_bumploc.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3denvar_bumploc_specified
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3denvar_bumploc_specified.yaml"
                       testoutput/3denvar_bumploc_specified.test
                  MPI     1
                  DEPENDS mpas_variational.x )

ecbuild_add_test( TARGET test_mpas_3dhybrid_bumpcov_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/3dhybrid_bumpcov_bumploc.yaml"
                       testoutput/3dhybrid_bumpcov_bumploc.test
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

ecbuild_add_test( TARGET test_mpas_4denvar_bumploc
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_variational.x testinput/4denvar_bumploc.yaml"
                       testoutput/4denvar_bumploc.test
                  MPI     1
                  DEPENDS mpas_variational.x )
