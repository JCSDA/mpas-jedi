/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "mpasjedi/ErrorCovarianceMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------

ErrorCovarianceMPAS::ErrorCovarianceMPAS(const GeometryMPAS & resol,
                                         const oops::Variables &,
                                         const eckit::Configuration & config,
                                         const StateMPAS &,
                                         const StateMPAS & ) {
  time_ = util::DateTime(config.getString("date"));
  mpas_b_setup_f90(keyErrCov_, config, resol.toFortran());
  oops::Log::trace() << "ErrorCovarianceMPAS created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovarianceMPAS::~ErrorCovarianceMPAS() {
  mpas_b_delete_f90(keyErrCov_);
  oops::Log::trace() << "ErrorCovarianceMPAS destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::linearize(const StateMPAS &,
                                    const GeometryMPAS & resol) {
  geom_.reset(new GeometryMPAS(resol));
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::multiply(const IncrementMPAS & dxin,
                                   IncrementMPAS & dxout) const {
  mpas_b_mult_f90(keyErrCov_, dxin.toFortran(),
                  dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::inverseMultiply(const IncrementMPAS & dxin,
                                           IncrementMPAS & dxout) const {
  mpas_b_invmult_f90(keyErrCov_, dxin.toFortran(),
                               dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::randomize(IncrementMPAS & dx) const {
  mpas_b_randomize_f90(keyErrCov_, dx.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::print(std::ostream & os) const {
  os << "ErrorCovarianceMPAS::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace mpas
