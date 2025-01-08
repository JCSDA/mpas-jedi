/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"
#include "mpasjedi/VariableChange/VariableChange.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------
  typedef VariableChangeParameters Parameters_;

VariableChange::VariableChange(const eckit::Configuration & config, const Geometry & geometry) {
  // Create the variable change
  VariableChangeParameters params;
  params.deserialize(config);
  variableChange_.reset(VariableChangeFactory::create(geometry,
                        params.variableChangeParametersWrapper.variableChangeParameters.value()));
}

// -------------------------------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars == x.variables()) {
    oops::Log::info() << "VariableChange::changeVar done (identity)" << std::endl;
    return;
  }

  oops::Log::trace() << "VariableChange::changeVar, vars" << vars << std::endl;
  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars == x.variables()) {
    oops::Log::info() << "VariableChange::changeVarInverse done (identity)" << std::endl;
    return;
  }

  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
