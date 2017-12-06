ecbuild_add_executable( TARGET  mpas_dfi.x
                        SOURCES mpasDFI.cc
                        LIBS    mpas
                       )

ecbuild_add_executable( TARGET  mpas_forecast.x
                        SOURCES mpasForecast.cc
                        LIBS    mpas
                       )

ecbuild_add_executable( TARGET  mpas_genpert.x
                        SOURCES mpasGenEnsPertB.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_4dvar.x
                        SOURCES mpas4DVar.cc
                        LIBS    mpas
                      )

ecbuild_add_executable( TARGET  mpas_makeobs.x
                        SOURCES mpasMakeObs.cc
                        LIBS    mpas
                      )

