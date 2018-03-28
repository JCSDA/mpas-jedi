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
//ORG  : keyConfig_(0), tstep_(), resol_(resol)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory"))
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
//BJJ: define setTrajectory for TlmIdMPAS .... from TlmMPAS
//ORG: void TlmIdMPAS::setTrajectory(const StateMPAS &, StateMPAS &, const ModelBiasMPAS &) {
void TlmIdMPAS::setTrajectory(const StateMPAS & xx, StateMPAS & xlr, const ModelBiasMPAS & bias) {
  oops::Log::trace() << "HERE IS TlmIdMPAS::setTrajectory" << std::endl;
// StateMPAS xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<double> zstat(15);
//  mpas_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "TlmIdMPAS trajectory at time " << xx.validTime() << std::endl;
//  for (unsigned int jj = 0; jj < 5; ++jj) {
//    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
//                       << ", RMS=" << zstat[3*jj+2] << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::initializeTL(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmIdMPAS::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::stepTL(IncrementMPAS & dx, const ModelBiasIncrementMPAS &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::finalizeTL(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmIdMPAS::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::initializeAD(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmIdMPAS::initializeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::stepAD(IncrementMPAS & dx, ModelBiasIncrementMPAS &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::finalizeAD(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmIdMPAS::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdMPAS::print(std::ostream & os) const {
  os << "MPAS IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
