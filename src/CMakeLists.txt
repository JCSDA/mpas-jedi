add_subdirectory( mpasjedi )
add_subdirectory( mains )

ecbuild_add_test( TARGET ${PROJECT_NAME}_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
