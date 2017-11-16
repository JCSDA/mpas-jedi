list( APPEND mpas_test_input
  testinput/interfaces.json
  testinput/test.json
  testinput/dirac.nicas.json
  compare.sh
)

list( APPEND qg_testoutput
  testoutput/dirac.nicas.test
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${mpas_test_input})
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

ecbuild_add_test( TARGET  test_mpas_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mpas )

ecbuild_add_test( TARGET test_mpas_dirac_nicas
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mpas_dirac.x testinput/dirac.nicas.json"
                       testoutput/dirac.nicas.test
                  DEPENDS mpas_dirac.x )

#ecbuild_add_test( TARGET  test_Interpolation
#//                  BOOST
##                  SOURCES executables/TestInterpolation.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    mpas )


