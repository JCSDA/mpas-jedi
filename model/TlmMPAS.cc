/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <vector>

#include "model/TlmMPAS.h"

#include "eckit/config/LocalConfiguration.h"
#include "Fortran.h"
#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "ModelMPAS.h"
#include "StateMPAS.h"
#include "MPASTraits.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/abor1_cpp.h"

namespace mpas {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<MPASTraits, TlmMPAS> makerMPASTLM_("MPASTLM");
// -----------------------------------------------------------------------------
TlmMPAS::TlmMPAS(const GeometryMPAS & resol,
                 const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory")),
    linvars_(std::vector<std::string>{"temperature", "pressure", "index_qv",
                              "uReconstructZonal", "uReconstructMeridional"})

{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  mpas_model_setup_f90(&configc, resol_.toFortran(), keyConfig_);

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
  mpas_model_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmMPAS::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::stepTL(IncrementMPAS & dx, const ModelBiasIncrementMPAS &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmMPAS: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("TlmMPAS: trajectory not available");
  }
  oops::Log::debug() << "TlmMPAS::stepTL fields in" << dx.fields() << std::endl;
  mpas_model_propagate_tl_f90(keyConfig_, dx.fields().toFortran(),
                              itra->second);
  oops::Log::debug() << "TlmMPAS::stepTL fields out" << dx.fields()
                     << std::endl;
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TlmMPAS::finalizeTL(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmMPAS::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::initializeAD(IncrementMPAS & dx) const {
  oops::Log::debug() << "TlmMPAS::initializeAD" << dx.fields() << std::endl;
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
  oops::Log::debug() << "TlmMPAS::stepAD fields in" << dx.fields() << std::endl;
  mpas_model_propagate_ad_f90(keyConfig_, dx.fields().toFortran(),
                              itra->second);
  oops::Log::debug() << "TlmMPAS::stepAD fields out" << dx.fields()
                     << std::endl;
}
// -----------------------------------------------------------------------------
void TlmMPAS::finalizeAD(IncrementMPAS & dx) const {
  mpas_model_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmMPAS::finalizeAD" << dx.fields() << std::endl;
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
