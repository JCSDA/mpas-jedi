
ecbuild_add_executable( TARGET  mpas_forecast.x
                        SOURCES mpasForecast.cc
                        LIBS    mpas
                       )

ecbuild_add_executable( TARGET  mpas_hofx.x
                        SOURCES mpasHofX.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_hofx3d.x
                        SOURCES mpasHofX3D.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_enshofx.x
                        SOURCES mpasEnsHofX.cc
                        LIBS    mpas
                      )

#ecbuild_add_executable( TARGET  mpas_4dvar.x
#                        SOURCES mpas4DVar.cc
#                        LIBS    mpas
#                      )

ecbuild_add_executable( TARGET  mpas_dirac.x
                        SOURCES mpasDirac.cc
                        LIBS	mpas
                      )

ecbuild_add_executable( TARGET  mpas_variational.x
                        SOURCES mpasVariational.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_parameters.x
                        SOURCES mpasEstimateParams.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_staticbinit.x
                        SOURCES mpasStaticBInit.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_gen_ens_pert_B.x
                        SOURCES mpasGenEnsPertB.cc
                        LIBS    mpas
                      )

