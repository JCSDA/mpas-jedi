list( APPEND oops_mpas_src_files
ErrorCovariance.h
Fields.cc
Fields.h
Fortran.h
Geometry.cc
Geometry.h
Increment.cc
Increment.h
State.cc
State.h
Traits.h
Variables.cc
Variables.h
random_f.cc
random_f.h
random_vectors_mod.F90
mpas_constants.f90
mpas_fields.F90
mpas_fields_interface.F90
mpas_geom_mod.F90
mpas_invent_state.f90
mpas_vars_mod.F90
)

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( cfdealloc C_F_NO_DEALLOC )
endif()

ecbuild_add_library( TARGET   mpas
                     DEFINITIONS ${cfdealloc}
                     SOURCES  ${oops_mpas_src_files}
                     LIBS     oops ${LAPACK_LIBRARIES} ${MPAS_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

