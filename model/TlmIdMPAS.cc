/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "TlmIdMPAS.h"

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "Fortran.h"
#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "MPASTraits.h"
#include "StateMPAS.h"
#include "util/DateTime.h"
#include "util/abor1_cpp.h"

namespace mpas {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<MPASTraits, TlmIdMPAS> makerMPASIdTLM_("MPASIdTLM");
// -----------------------------------------------------------------------------
TlmIdMPAS::TlmIdMPAS(const GeometryMPAS & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol)
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  mpas_model_setup_f90(&configc, resol_.toFortran(), keyConfig_);

  oops::Log::trace() << "TlmIdMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmIdMPAS::~TlmIdMPAS() {
  mpas_model_delete_f90(keyConfig_);
  oops::Log::trace() << "TlmIdMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::setTrajectory(const StateMPAS &, StateMPAS &, const ModelBiasMPAS &) {}
// -----------------------------------------------------------------------------
void TlmIdMPAS::initializeTL(IncrementMPAS & dx) const {
  dx.activateModel();
  mpas_model_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmIdMPAS::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::stepTL(IncrementMPAS & dx, const ModelBiasIncrementMPAS &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::finalizeTL(IncrementMPAS & dx) const {
  dx.deactivateModel();
  oops::Log::debug() << "TlmIdMPAS::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::initializeAD(IncrementMPAS & dx) const {
  dx.activateModel();
  oops::Log::debug() << "TlmIdMPAS::initializeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::stepAD(IncrementMPAS & dx, ModelBiasIncrementMPAS &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::finalizeAD(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  dx.deactivateModel();
  oops::Log::debug() << "TlmIdMPAS::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::print(std::ostream & os) const {
  os << "MPAS IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
