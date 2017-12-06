ecbuild_add_executable( TARGET  lfric_dfi.x
                        SOURCES lfricDFI.cc
                        LIBS    lfric
                       )

ecbuild_add_executable( TARGET  lfric_forecast.x
                        SOURCES lfricForecast.cc
                        LIBS    lfric
                       )

ecbuild_add_executable( TARGET  lfric_genpert.x
                        SOURCES lfricGenEnsPertB.cc
                        LIBS    lfric
                      )

ecbuild_add_executable( TARGET  lfric_4dvar.x
                        SOURCES lfric4DVar.cc
                        LIBS    lfric
                      )

ecbuild_add_executable( TARGET  lfric_makeobs.x
                        SOURCES lfricMakeObs.cc
                        LIBS    lfric
                      )

