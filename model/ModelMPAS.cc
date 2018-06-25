/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/ModelMPAS.h"

#include "oops/util/Logger.h"
#include "ModelBiasMPAS.h"
#include "FieldsMPAS.h"
#include "Fortran.h"
#include "GeometryMPAS.h"
#include "StateMPAS.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"

namespace mpas {
// -----------------------------------------------------------------------------
ModelMPAS::ModelMPAS(const GeometryMPAS & resol,
                     const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol)
{
  oops::Log::trace() << "ModelMPAS::ModelMPAS" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  oops::Log::trace() << "ModelMPAS::tstep_" << tstep_ << std::endl;
  const eckit::Configuration * configc = &model;
  mpas_model_setup_f90(&configc, geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelMPAS::~ModelMPAS() {
  mpas_model_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::initialize(StateMPAS & xx) const {
  mpas_model_prepare_integration_f90(keyConfig_, xx.fields().toFortran());
  oops::Log::debug() << "ModelMPAS::initialize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::step(StateMPAS & xx, const ModelBiasMPAS &) const {
  oops::Log::debug() << "ModelMPAS::step fields in" << xx.fields() << std::endl;
  mpas_model_propagate_f90(keyConfig_, xx.fields().toFortran());
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelMPAS::step fields out" << xx.fields()
                     << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::finalize(StateMPAS & xx) const {
  oops::Log::debug() << "ModelMPAS::finalize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
int ModelMPAS::saveTrajectory(StateMPAS & xx, const ModelBiasMPAS &) const {
  int ftraj = 0;
  oops::Log::debug() << "ModelMPAS::saveTrajectory fields in" << xx.fields()
                     << std::endl;
  mpas_model_prop_traj_f90(keyConfig_, xx.fields().toFortran(), ftraj);
  ASSERT(ftraj != 0);
  oops::Log::debug() << "ModelMPAS::saveTrajectory fields out" << xx.fields()
                     << std::endl;
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelMPAS::print(std::ostream & os) const {
  os << "ModelMPAS::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
