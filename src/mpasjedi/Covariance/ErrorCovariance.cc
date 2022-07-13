/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Covariance/ErrorCovariance.h"
#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/State/State.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------

ErrorCovariance::ErrorCovariance(const Geometry & resol,
                                         const oops::Variables &,
                                         const eckit::Configuration & config,
                                         const State &,
                                         const State & ) {
  time_ = util::DateTime(config.getString("date"));
  mpas_b_setup_f90(keyErrCov_, config, resol.toFortran());
  oops::Log::trace() << "ErrorCovariance created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovariance::~ErrorCovariance() {
  mpas_b_delete_f90(keyErrCov_);
  oops::Log::trace() << "ErrorCovariance destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovariance::linearize(const State &,
                                const Geometry & resol)  {
  ABORT("ErrorCovariance::linearize not implemented");
}

// -----------------------------------------------------------------------------

void ErrorCovariance::multiply(const Increment & dxin,
                                   Increment & dxout) const {
  mpas_b_mult_f90(keyErrCov_, dxin.toFortran(),
                  dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::inverseMultiply(const Increment & dxin,
                                           Increment & dxout) const {
  mpas_b_invmult_f90(keyErrCov_, dxin.toFortran(),
                               dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::randomize(Increment & dx) const {
  mpas_b_randomize_f90(keyErrCov_, dx.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::print(std::ostream & os) const {
  os << "ErrorCovariance::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace mpas
