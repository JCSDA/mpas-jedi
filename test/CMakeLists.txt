list( APPEND mpas_test_input
  testinput/3dvar.json
  testinput/3dfgat.json
  testinput/4densvar.json
  testinput/4dvar.alpha.json
  testinput/4dforcing.json
  testinput/4dhybrid.json
  testinput/4dsaddlepoint.json
  testinput/4dvar.dripcg.json
  testinput/4dvar.drpcglmp.json
  testinput/4dvar.drpfom.json
  testinput/4dvar.drplanczos.json
  testinput/4dvar.ipcg.json
  testinput/4dvar.obsbias.json
  testinput/4dvar.rpcg.json
  testinput/dfi.json
  testinput/forecast.json
  testinput/genenspert.json
  testinput/interfaces.json
  testinput/makeobs3d.json
  testinput/makeobs4d.json
  testinput/makeobs4d12h.json
  testinput/makeobsbias.json
  testinput/test.json
  testinput/truth.json
  testinput/test_op_obs.json
  compare.sh
  ../statics/CAM_ABS_DATA.DBL
  ../statics/CAM_AEROPT_DATA.DBL
  ../statics/GENPARM.TBL
  ../statics/LANDUSE.TBL
  ../statics/namelist.atmosphere
  ../statics/OZONE_DAT.TBL
  ../statics/OZONE_LAT.TBL
  ../statics/OZONE_PLEV.TBL
  ../statics/RRTMG_LW_DATA.DBL
  ../statics/RRTMG_SW_DATA.DBL
  ../statics/SOILPARM.TBL
  ../statics/stream_list.atmosphere.diagnostics
  ../statics/stream_list.atmosphere.output
  ../statics/stream_list.atmosphere.surface
  ../statics/VEGPARM.TBL
  ../statics/streams.atmosphere
)

list( APPEND mpas_testoutput
  testoutput/3dvar.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4dvar.alpha.test
  testoutput/4dforcing.test
  testoutput/4dhybrid.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar.dripcg.test
  testoutput/4dvar.drpcglmp.test
  testoutput/4dvar.drpfom.test
  testoutput/4dvar.drplanczos.test
  testoutput/4dvar.ipcg.test
  testoutput/4dvar.obsbias.test
  testoutput/4dvar.rpcg.test
  testoutput/dfi.test
  testoutput/forecast.test
  testoutput/genenspert.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/truth.test
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

#ecbuild_add_test( TARGET test_mpas_truth
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_forecast.x testinput/truth.json"
#                       testoutput/truth.test
#                  DEPENDS mpas_forecast.x )
#
#ecbuild_add_test( TARGET test_mpas_forecast
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_forecast.x testinput/forecast.json"
#                       testoutput/forecast.test
#                  DEPENDS mpas_forecast.x )
#
#ecbuild_add_test( TARGET test_mpas_makeobs3d
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_makeobs.x testinput/makeobs3d.json"
#                       testoutput/makeobs3d.test
#                  DEPENDS mpas_makeobs.x )
#
#ecbuild_add_test( TARGET test_mpas_makeobs4d
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_makeobs.x testinput/makeobs4d.json"
#                       testoutput/makeobs4d.test
#                  DEPENDS mpas_makeobs.x )
#
#ecbuild_add_test( TARGET test_mpas_makeobs4d12h
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_makeobs.x testinput/makeobs4d12h.json"
#                       testoutput/makeobs4d12h.test
#                  DEPENDS mpas_makeobs.x )
#
#ecbuild_add_test( TARGET test_mpas_makeobsbias
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_makeobs.x testinput/makeobsbias.json"
#                       testoutput/makeobsbias.test
#                  DEPENDS mpas_makeobs.x )

# Test interface classes with MPAS

ecbuild_add_test( TARGET  test_mpas_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mpas )

ecbuild_add_test( TARGET  test_mpas_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mpas )
