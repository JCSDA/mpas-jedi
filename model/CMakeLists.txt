list( APPEND oops_mpas_src_files
ErrorCovarianceMPAS.cc
ErrorCovarianceMPAS.h
VariablesMPAS.cc
VariablesMPAS.h
FieldsMPAS.cc
FieldsMPAS.h
Fortran.h
GeometryMPAS.cc
GeometryMPAS.h
IncrementMPAS.cc
IncrementMPAS.h
MPASTraits.h
LocalizationMatrixMPAS.cc
LocalizationMatrixMPAS.h
ModelBiasCovarianceMPAS.h
ModelBiasIncrementMPAS.h
ModelBiasMPAS.h
ModelMPAS.cc
ModelMPAS.h
RunMPAS.cc
RunMPAS.h
StateMPAS.cc
StateMPAS.h
TlmIdMPAS.cc
TlmIdMPAS.h
TlmMPAS.cc
TlmMPAS.h
#VariablesMPAS.cc
#VariablesMPAS.h
mpas_covariance_interface.F90
mpas_covariance_mod.F90
mpas_fields_interface.F90
mpas_fields_mod.F90
mpas_geom_interface.F90
mpas_geom_mod.F90
mpas_kinds_mod.F90
mpas_model_interface.F90
mpas_model_mod.F90
mpas_trajectories.F90
tools_nc.f90
mpas_vars_mod.F90
mpas2da_mod.F90
mpas_run_interface.F90
#mpas_fields_interp_mod.F90
)

set ( MPAS_LIBRARIES $ENV{MPAS_LIBRARIES} )
set ( MPAS_INCLUDES $ENV{MPAS_INCLUDES} )

ecbuild_add_library( TARGET   mpas
                     SOURCES  ${oops_mpas_src_files}
                     LIBS     oops ufo ${MPAS_LIBRARIES}
                     INCLUDES ${MPAS_INCLUDES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

