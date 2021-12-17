/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/StateMPAS.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const GeometryMPAS & geom, const Parameters_ & params)
  : geom_(new GeometryMPAS(geom)), params_(params), linearVariableChange_() {}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::setTrajectory(const StateMPAS & xbg, const StateMPAS & xfg) {
  oops::Log::trace() << "LinearVariableChange::setTrajectory starting" << std::endl;
  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xbg, xfg, *geom_,
             params_.linearVariableChangeParametersWrapper.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiply(IncrementMPAS & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiply starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars <= dx.variables()) {
    dx.updateVars(vars);
    oops::Log::trace() << "LinearVariableChange::multiply done (identity)" << std::endl;
    return;
  }

  // Create output increment
  IncrementMPAS dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiply(dx, dxout);

  // Update variable names
  dx.updateVars(vars);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(IncrementMPAS & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverse starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars <= dx.variables()) {
    dx.updateVars(vars);
    oops::Log::trace() << "LinearVariableChange::multiplyInverse done (identity)" << std::endl;
    return;
  }

  // Create output increment
  IncrementMPAS dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverse(dx, dxout);

  // Update variable names
  dx.updateVars(vars);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(IncrementMPAS & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars <= dx.variables()) {
    dx.updateVars(vars);
    oops::Log::trace() << "LinearVariableChange::multiplyAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  IncrementMPAS dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyAD(dx, dxout);

  // Update variable names
  dx.updateVars(vars);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(IncrementMPAS & dx,
                                             const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars <= dx.variables()) {
    dx.updateVars(vars);
    oops::Log::trace() << "LinearVariableChange::multiplyInverseAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  IncrementMPAS dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverseAD(dx, dxout);

  // Update variable names
  dx.updateVars(vars);

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
