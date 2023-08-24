/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {
  class Variables;
}

namespace mpas {
// -----------------------------------------------------------------------------
/// Model Parameters Class.

class ModelParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelParameters, Parameters)

 public:
/// The property 'name' is already part of the default schema
/// Model State variables to process
  oops::RequiredParameter<oops::Variables> vars{"model variables", this};
/// Model time step
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
};
// -----------------------------------------------------------------------------

}  // namespace mpas
