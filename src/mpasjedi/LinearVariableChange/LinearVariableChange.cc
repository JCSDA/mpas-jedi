/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom, const Parameters_ & params)
  : geom_(new Geometry(geom)), params_(params), linearVariableChange_() {}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::setTrajectory(const State & xbg, const State & xfg) {
  oops::Log::trace() << "LinearVariableChange::setTrajectory starting" << std::endl;
  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xbg, xfg, *geom_,
             params_.linearVariableChangeParametersWrapper.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiply(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiply starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::multiply done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiply(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverse starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::multiplyInverse done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverse(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::multiplyAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyAD(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(Increment & dx,
                                             const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::multiplyInverseAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverseAD(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << classname() << " linear variable change";
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
