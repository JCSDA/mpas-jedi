/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/interface/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "mpasjedi/Model/ModelParameters.h"

namespace oops {
  class Variables;
}

namespace mpas {
// -----------------------------------------------------------------------------
/// Tlm Parameters Class.

class TlmParameters : public oops::LinearModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(TlmParameters, LinearModelParametersBase)

 public:
/// Tlm Model State variables to process
  oops::RequiredParameter<oops::Variables> tlmvars{"tlm variables", this};
/// Tlm Model time step
  oops::RequiredParameter<util::Duration> tlmtstep{"tstep", this};
  ModelParameters tlmparams {this};  // includes model parameters
};
// -----------------------------------------------------------------------------

}  // namespace mpas
