/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Fortran.h"
#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/ModelBiasMPAS.h"
#include "mpasjedi/ModelMPAS.h"
#include "mpasjedi/StateMPAS.h"

namespace mpas {
// -----------------------------------------------------------------------------
static oops::interface::ModelMaker<MPASTraits, ModelMPAS> makermodel_("MPAS");
// -----------------------------------------------------------------------------
ModelMPAS::ModelMPAS(const GeometryMPAS & resol,
                     const eckit::Configuration & model)
  : keyModel_(0), tstep_(0), geom_(resol),
    vars_(model, "model variables")
{
  oops::Log::trace() << "ModelMPAS::ModelMPAS" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  oops::Log::trace() << "ModelMPAS::tstep_" << tstep_ << std::endl;
  mpas_model_setup_f90(model, geom_.toFortran(), keyModel_);
  oops::Log::trace() << "ModelMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelMPAS::~ModelMPAS() {
  mpas_model_delete_f90(keyModel_);
  oops::Log::trace() << "ModelMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::initialize(StateMPAS & xx) const {
  mpas_model_prepare_integration_f90(keyModel_, xx.toFortran());
  oops::Log::debug() << "ModelMPAS::initialize" << xx << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::step(StateMPAS & xx, const ModelBiasMPAS &) const {
  oops::Log::debug() << "ModelMPAS::step state in" << xx << std::endl;
  mpas_model_propagate_f90(keyModel_, xx.toFortran());
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelMPAS::step state out" << xx << std::endl;
}
// -----------------------------------------------------------------------------
void ModelMPAS::finalize(StateMPAS & xx) const {
  oops::Log::debug() << "ModelMPAS::finalize" << xx << std::endl;
}
// -----------------------------------------------------------------------------
int ModelMPAS::saveTrajectory(StateMPAS & xx, const ModelBiasMPAS &) const {
  int ftraj = 0;
  oops::Log::debug() << "ModelMPAS::saveTrajectory state in" << xx << std::endl;
  mpas_model_prop_traj_f90(keyModel_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  oops::Log::debug() << "ModelMPAS::saveTrajectory state out" << xx <<std::endl;
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelMPAS::print(std::ostream & os) const {
  os << "ModelMPAS::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
