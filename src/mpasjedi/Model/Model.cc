/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Model/Model.h"
#include "mpasjedi/Model/ModelParameters.h"
#include "mpasjedi/ModelBias/ModelBias.h"
#include "mpasjedi/State/State.h"

namespace mpas {
// -----------------------------------------------------------------------------
static oops::interface::ModelMaker<Traits, Model> makermodel_("MPAS");
// -----------------------------------------------------------------------------
Model::Model(const Geometry & resol,
                     const eckit::Configuration & config)
  : keyModel_(0), tstep_(0), geom_(resol), vars_(config, "model variables")
{
  oops::Log::trace() << "Model::Model" << std::endl;
  ModelParameters params;
  params.deserialize(config);
  tstep_ = util::Duration(config.getString("tstep"));
  oops::Log::trace() << "Model::tstep_" << tstep_ << std::endl;
  mpas_model_setup_f90(params.toConfiguration(), geom_.toFortran(), keyModel_);
  oops::Log::trace() << "Model created" << std::endl;
}
// -----------------------------------------------------------------------------
Model::~Model() {
  mpas_model_delete_f90(keyModel_);
  oops::Log::trace() << "Model destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void Model::initialize(State & xx) const {
  mpas_model_prepare_integration_f90(keyModel_, xx.toFortran());
  oops::Log::debug() << "Model::initialize" << xx << std::endl;
}
// -----------------------------------------------------------------------------
void Model::step(State & xx, const ModelBias &) const {
  oops::Log::debug() << "Model::step state in" << xx << std::endl;
  mpas_model_propagate_f90(keyModel_, xx.toFortran());
  xx.validTime() += tstep_;
  oops::Log::debug() << "Model::step state out" << xx << std::endl;
}
// -----------------------------------------------------------------------------
void Model::finalize(State & xx) const {
  oops::Log::debug() << "Model::finalize" << xx << std::endl;
}
// -----------------------------------------------------------------------------
int Model::saveTrajectory(State & xx, const ModelBias &) const {
  int ftraj = 0;
  oops::Log::debug() << "Model::saveTrajectory state in" << xx << std::endl;
  mpas_model_prop_traj_f90(keyModel_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  oops::Log::debug() << "Model::saveTrajectory state out" << xx <<std::endl;
  return ftraj;
}
// -----------------------------------------------------------------------------
void Model::print(std::ostream & os) const {
  os << "Model::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
