/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"
#include "mpasjedi/TlmMPAS.h"

namespace mpas {

// -----------------------------------------------------------------------------
static oops::interface::LinearModelMaker<MPASTraits, TlmMPAS> makerMPASTLM_("MPASTLM");
// -----------------------------------------------------------------------------
TlmMPAS::TlmMPAS(const GeometryMPAS & resol,
                 const Parameters_ & params)
  : keyConfig_(0), tstep_(params.tlmtstep), resol_(resol), traj_(),
    lrmodel_(resol_, params.tlmparams),
    linvars_(params.tlmvars)

{
  tstep_ = util::Duration(params.tlmtstep);

  mpas_model_setup_f90(params.toConfiguration(), resol_.toFortran(), keyConfig_);

  oops::Log::trace() << "TlmMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmMPAS::~TlmMPAS() {
  mpas_model_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    mpas_model_wipe_traj_f90(jtra->second);
  }
  oops::Log::trace() << "TlmMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::setTrajectory(const StateMPAS & xx, StateMPAS & xlr,
                            const ModelBiasMPAS & bias) {
// StateMPAS xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<double> zstat(15);
//  mpas_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "TlmMPAS trajectory at time " << xx.validTime()
                     << std::endl;
  for (unsigned int jj = 0; jj < 5; ++jj) {
    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
                       << ", RMS=" << zstat[3*jj+2] << std::endl;
  }
// should be in print method
}
// -----------------------------------------------------------------------------
void TlmMPAS::initializeTL(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_tl_f90(keyConfig_, dx.toFortran());
  oops::Log::debug() << "TlmMPAS::initializeTL" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::stepTL(IncrementMPAS & dx, const ModelBiasIncrementMPAS &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmMPAS: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("TlmMPAS: trajectory not available");
  }
  oops::Log::debug() << "TlmMPAS::stepTL increment in" << dx << std::endl;
  mpas_model_propagate_tl_f90(keyConfig_, dx.toFortran(),
                              itra->second);
  oops::Log::debug() << "TlmMPAS::stepTL increment out"<< dx << std::endl;
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TlmMPAS::finalizeTL(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmMPAS::finalizeTL" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::initializeAD(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmMPAS::initializeAD" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::stepAD(IncrementMPAS & dx, ModelBiasIncrementMPAS &) const {
  dx.validTime() -= tstep_;
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmMPAS: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("TlmMPAS: trajectory not available");
  }
  oops::Log::debug() << "TlmMPAS::stepAD increment in" << dx << std::endl;
  mpas_model_propagate_ad_f90(keyConfig_, dx.toFortran(),
                              itra->second);
  oops::Log::debug() << "TlmMPAS::stepAD increment out" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::finalizeAD(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_ad_f90(keyConfig_, dx.toFortran());
  oops::Log::debug() << "TlmMPAS::finalizeAD" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::print(std::ostream & os) const {
  os << "MPAS TLM Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "MPAS TLM Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace mpas
