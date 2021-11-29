
ecbuild_add_executable( TARGET  ${PROJECT_NAME}_dirac.x
                        SOURCES mpasDirac.cc
                        LIBS    ${PROJECT_NAME} saber )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_enshofx.x
                        SOURCES mpasEnsHofX.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_eda.x
                        SOURCES mpasEDA.cc
                        LIBS    ${PROJECT_NAME} saber )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_forecast.x
                        SOURCES mpasForecast.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_gen_ens_pert_B.x
                        SOURCES mpasGenEnsPertB.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_hofx.x
                        SOURCES mpasHofX.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_hofx3d.x
                        SOURCES mpasHofX3D.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_error_covariance_training.x
                        SOURCES mpasErrorCovarianceTraining.cc
                        LIBS    ${PROJECT_NAME} saber )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_rtpp.x
                        SOURCES mpasRTPP.cc
                        LIBS	${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_staticbinit.x
                        SOURCES mpasStaticBInit.cc
                        LIBS    ${PROJECT_NAME} )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_variational.x
                        SOURCES mpasVariational.cc
                        LIBS    ${PROJECT_NAME} saber )

ecbuild_add_executable( TARGET  ${PROJECT_NAME}_convertstate.x
                        SOURCES mpasConvertState.cc
                        LIBS    ${PROJECT_NAME} saber )
