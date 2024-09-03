list( APPEND mpasjedi_src_files
    Covariance/ErrorCovariance.cc
    Covariance/ErrorCovariance.h
    Covariance/ErrorCovariance.interface.h
    Covariance/mpas_covariance_interface.F90
    Covariance/mpas_covariance_mod.F90
    Fields/mpas_fields_mod.F90
    Fortran.h
    Geometry/Geometry.cc
    Geometry/Geometry.h
    Geometry/GeometryParameters.h
    Geometry/Geometry.interface.h
    Geometry/mpas_geom_interface.F90
    Geometry/mpas_geom_mod.F90
    GeometryIterator/GeometryIterator.cc
    GeometryIterator/GeometryIterator.h
    GeometryIterator/GeometryIterator.interface.h
    GeometryIterator/iter.F90
    GeometryIterator/iter2d.F90
    GeometryIterator/iter3d.F90
    GeometryIterator/mpas_geom_iter_mod.F90
    GeometryIterator/mpas_geom_iter_interface.F90
    Increment/Increment.cc
    Increment/Increment.h
    Increment/IncrementParameters.h
    Increment/Increment.interface.h
    Increment/mpas_increment_interface_mod.F90
    Increment/mpas_increment_mod.F90
    LinearVariableChange/LinearVariableChange.cc
    LinearVariableChange/LinearVariableChange.h
    LinearVariableChange/Base/LinearVariableChangeBase.cc
    LinearVariableChange/Base/LinearVariableChangeBase.h
    LinearVariableChange/Control2Analysis/mpasjedi_linvarcha_c2a_interface.F90
    LinearVariableChange/Control2Analysis/mpasjedi_linvarcha_c2a_mod.F90
    LinearVariableChange/Control2Analysis/LinVarChaC2A.cc
    LinearVariableChange/Control2Analysis/LinVarChaC2A.h
    LinearVariableChange/Control2Analysis/LinVarChaC2A.interface.h
    LinearVariableChange/Model2GeoVars/mpasjedi_lvc_model2geovars_interface.F90
    LinearVariableChange/Model2GeoVars/mpasjedi_lvc_model2geovars_mod.F90
    LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.cc
    LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.h
    LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.interface.h
    Model/Model.cc
    Model/Model.h
    Model/ModelParameters.h
    Model/Model.interface.h
    Model/mpas_model_interface.F90
    Model/mpas_model_mod.F90
    Model/mpas_trajectories.F90
    ModelBias/ModelBiasCovariance.h
    ModelBias/ModelBiasIncrement.h
    ModelBias/ModelBias.h
    ModelData/ModelData.h
    Saca/mpas_saca_interface_mod.F90
    Saca/module_mp_thompson_cldfra3.F90
    State/State.cc
    State/State.h
    State/State.interface.h
    State/StateParameters.h 
    State/mpas_state_interface_mod.F90
    State/mpas_state_mod.F90
    Tlm/Tlm.cc
    Tlm/Tlm.h
    Traits.h
    VariableChange/VariableChange.cc
    VariableChange/VariableChange.h
    VariableChange/Base/VariableChangeBase.cc
    VariableChange/Base/VariableChangeBase.h
    VariableChange/Model2GeoVars/mpasjedi_vc_model2geovars_interface.F90
    VariableChange/Model2GeoVars/mpasjedi_vc_model2geovars_mod.F90
    VariableChange/Model2GeoVars/VarChaModel2GeoVars.cc
    VariableChange/Model2GeoVars/VarChaModel2GeoVars.h
    VariableChange/Model2GeoVars/VarChaModel2GeoVars.interface.h
    Utilities/mpas2ufo_vars_mod.F90
    Utilities/mpas4da_mod.F90
    Utilities/mpas_constants_mod.F90
    Utilities/mpas_kinds_mod.F90
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
