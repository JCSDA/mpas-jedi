list( APPEND mpasjedi_src_files
    ErrorCovarianceMPAS.cc
    ErrorCovarianceMPAS.h
    Fortran.h
    GeometryMPAS.cc
    GeometryMPAS.h
    IncrementMPAS.cc
    IncrementMPAS.h
    IncrementMPASFortran.h
    MPASTraits.h
    ModelBiasCovarianceMPAS.h
    ModelBiasIncrementMPAS.h
    ModelBiasMPAS.h
    ModelMPAS.cc
    ModelMPAS.h
    StateMPAS.cc
    StateMPAS.h
    StateMPASFortran.h
    TlmMPAS.cc
    TlmMPAS.h
    mpas_constants_mod.F90
    mpas_covariance_interface.F90
    mpas_covariance_mod.F90
    mpas_fields_mod.F90
    mpas_geom_interface.F90
    mpas_geom_mod.F90
    mpas_increment_interface_mod.F90
    mpas_increment_mod.F90
    mpas_model_interface.F90
    mpas_model_mod.F90
    mpas_state_interface_mod.F90
    mpas_state_mod.F90
    mpas_trajectories.F90
    mpas2ufo_vars_mod.F90
    mpas4da_mod.F90
    mpas_kinds_mod.F90
    getvalues/mpasjedi_getvalues_mod.F90
    getvalues/mpasjedi_lineargetvalues_mod.F90
    getvalues/GetValues.cc
    getvalues/GetValues.h
    getvalues/LinearGetValues.cc
    getvalues/LinearGetValues.h
    getvalues/GetValues.interface.h
    getvalues/GetValues.interface.F90
    getvalues/LinearGetValues.interface.h
    getvalues/LinearGetValues.interface.F90
    getvalues/mpasjedi_unstructured_interp_mod.F90
    VariableChanges/Control2Analysis/mpasjedi_linvarcha_c2a_interface.F90
    VariableChanges/Control2Analysis/mpasjedi_linvarcha_c2a_mod.F90
    VariableChanges/Control2Analysis/LinVarChaC2A.cc
    VariableChanges/Control2Analysis/LinVarChaC2A.h
    VariableChanges/Control2Analysis/LinVarChaC2A.interface.h
    VariableChanges/Model2GeoVars/mpasjedi_lvc_model2geovars_interface.F90
    VariableChanges/Model2GeoVars/mpasjedi_lvc_model2geovars_mod.F90
    VariableChanges/Model2GeoVars/mpasjedi_vc_model2geovars_interface.F90
    VariableChanges/Model2GeoVars/mpasjedi_vc_model2geovars_mod.F90
    VariableChanges/Model2GeoVars/LinVarChaModel2GeoVars.cc
    VariableChanges/Model2GeoVars/LinVarChaModel2GeoVars.h
    VariableChanges/Model2GeoVars/LinVarChaModel2GeoVars.interface.h
    VariableChanges/Model2GeoVars/VarChaModel2GeoVars.cc
    VariableChanges/Model2GeoVars/VarChaModel2GeoVars.h
    VariableChanges/Model2GeoVars/VarChaModel2GeoVars.interface.h
    )

#Create target
include(GNUInstallDirs)
ecbuild_add_library( TARGET ${PROJECT_NAME}
                     SOURCES ${mpasjedi_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
                     )

#Link dependencies
if(OpenMP_FOUND)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( ${PROJECT_NAME} PUBLIC Boost::boost )
target_link_libraries( ${PROJECT_NAME} PUBLIC MPAS::core::atmosphere )
target_link_libraries( ${PROJECT_NAME} PUBLIC oops )
target_link_libraries( ${PROJECT_NAME} PUBLIC saber )
target_link_libraries( ${PROJECT_NAME} PUBLIC ioda )
target_link_libraries( ${PROJECT_NAME} PUBLIC ufo )

#Configure include directory layout for build-tree to match install-tree
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

#Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                                $<INSTALL_INTERFACE:${MODULE_DIR}>)
