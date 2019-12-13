list( APPEND oops_mpas_src_files
ErrorCovarianceMPAS.cc
ErrorCovarianceMPAS.h
Fortran.h
GeometryMPAS.cc
GeometryMPAS.h
GetValuesTrajMPAS.cc
GetValuesTrajMPAS.h
IncrementMPAS.cc
IncrementMPAS.h
IncrementMPASFortran.h
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
StateMPASFortran.h
TlmIdMPAS.cc
TlmIdMPAS.h
TlmMPAS.cc
TlmMPAS.h
VarChangeMPAS.cc
VarChangeMPAS.h
mpas_constants_mod.F90
mpas_covariance_interface.F90
mpas_covariance_mod.F90
mpas_field_utils_mod.F90
mpas_geom_interface.F90
mpas_geom_mod.F90
mpas_getvaltraj_mod.F90
mpas_increment_interface_mod.F90
mpas_increment_mod.F90
mpas_model_interface.F90
mpas_model_mod.F90
mpas_run_mod.F90
mpas_state_interface_mod.F90
mpas_state_mod.F90
mpas_trajectories.F90
mpas_varchange_interface.F90
mpas_varchange_mod.F90
mpas2ufo_vars_mod.F90
mpas4da_mod.F90
mpas_kinds_mod.F90
instantiateLocalizationFactory.h
instantiateMPASVarChangeFactory.h
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

