/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/ErrorCovarianceMPAS.h"

#include <cmath>

#include "oops/util/Logger.h"
#include "Fortran.h"
#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "StateMPAS.h"
#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------

ErrorCovarianceMPAS::ErrorCovarianceMPAS(const GeometryMPAS & resol,
                                         const oops::Variables &,
                                         const eckit::Configuration & conf,
                                         const StateMPAS &,
                                         const StateMPAS & ) {
  time_ = util::DateTime(conf.getString("date"));
  const eckit::Configuration * configc = &conf;
  mpas_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
  oops::Log::trace() << "ErrorCovarianceMPAS created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovarianceMPAS::~ErrorCovarianceMPAS() {
  mpas_b_delete_f90(keyFtnConfig_);
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
  mpas_b_mult_f90(keyFtnConfig_, dxin.toFortran(),
                  dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::inverseMultiply(const IncrementMPAS & dxin,
                                           IncrementMPAS & dxout) const {
  mpas_b_invmult_f90(keyFtnConfig_, dxin.toFortran(),
                               dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::randomize(IncrementMPAS & dx) const {
  mpas_b_randomize_f90(keyFtnConfig_, dx.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceMPAS::print(std::ostream & os) const {
  os << "ErrorCovarianceMPAS::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace mpas
