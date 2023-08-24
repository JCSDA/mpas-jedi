/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const eckit::Configuration & config)
  : geom_(geom), linearVariableChange_() {
  params_.deserialize(config);
}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg, const oops::Variables & vars) {
  oops::Log::trace() << "LinearVariableChange::changeVarTraj starting" << std::endl;
  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xfg, xfg, geom_,
             params_.linearVariableChangeParametersWrapper.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::changeVarTraj done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarTL starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarTL done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarTL(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarTL done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarInverseTL(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarAD(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                             const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarInverseAD(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << classname() << " linear variable change";
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
