/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/Model/Model.interface.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Tlm/Tlm.h"

namespace mpas {

// -----------------------------------------------------------------------------
static oops::interface::LinearModelMaker<Traits, Tlm> makerMPASTLM_("MPASTLM");
// -----------------------------------------------------------------------------
Tlm::Tlm(const Geometry & resol,
         const eckit::Configuration & config)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, config)
{
  tstep_ = util::Duration(config.getString("tstep"));

  mpas_model_setup_f90(config, resol_.toFortran(), keyConfig_);

  oops::Log::trace() << "Tlm created" << std::endl;
}
// -----------------------------------------------------------------------------
Tlm::~Tlm() {
  mpas_model_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    mpas_model_wipe_traj_f90(jtra->second);
  }
  oops::Log::trace() << "Tlm destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::setTrajectory(const State & xx, State & xlr,
                            const ModelBias & bias) {
// State xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<real_type> zstat(15);
//  mpas_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "Tlm trajectory at time " << xx.validTime()
                     << std::endl;
  for (unsigned int jj = 0; jj < 5; ++jj) {
    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
                       << ", RMS=" << zstat[3*jj+2] << std::endl;
  }
// should be in print method
}
// -----------------------------------------------------------------------------
void Tlm::initializeTL(Increment & dx) const {
  mpas_model_prepare_integration_tl_f90(keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::initializeTL" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::stepTL(Increment & dx, const ModelBiasIncrement &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }
  oops::Log::debug() << "Tlm::stepTL increment in" << dx << std::endl;
  mpas_model_propagate_tl_f90(keyConfig_, dx.toFortran(),
                              itra->second);
  oops::Log::debug() << "Tlm::stepTL increment out"<< dx << std::endl;
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void Tlm::finalizeTL(Increment & dx) const {
  oops::Log::debug() << "Tlm::finalizeTL" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::initializeAD(Increment & dx) const {
  oops::Log::debug() << "Tlm::initializeAD" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::stepAD(Increment & dx, ModelBiasIncrement &) const {
  dx.validTime() -= tstep_;
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }
  oops::Log::debug() << "Tlm::stepAD increment in" << dx << std::endl;
  mpas_model_propagate_ad_f90(keyConfig_, dx.toFortran(),
                              itra->second);
  oops::Log::debug() << "Tlm::stepAD increment out" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::finalizeAD(Increment & dx) const {
  mpas_model_prepare_integration_ad_f90(keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::finalizeAD" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::print(std::ostream & os) const {
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
