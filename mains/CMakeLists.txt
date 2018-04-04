
ecbuild_add_executable( TARGET  mpas_forecast.x
                        SOURCES mpasForecast.cc
                        LIBS    mpas
                       )

ecbuild_add_executable( TARGET  mpas_hofx.x
                        SOURCES mpasHofX.cc
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

