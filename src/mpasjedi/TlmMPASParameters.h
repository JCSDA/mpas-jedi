/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPASJEDI_TLMMPASPARAMETERS_H_
#define MPASJEDI_TLMMPASPARAMETERS_H_

#include "oops/interface/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "mpasjedi/ModelMPASParameters.h"

namespace oops {
  class Variables;
}

namespace mpas {
// -----------------------------------------------------------------------------
/// Tlm Parameters Class.

class TlmMPASParameters : public oops::LinearModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(TlmMPASParameters, LinearModelParametersBase)

 public:
/// Tlm Model State variables to process
  oops::RequiredParameter<oops::Variables> tlmvars{"tlm variables", this};
/// Tlm Model time step
  oops::RequiredParameter<util::Duration> tlmtstep{"tstep", this};
  ModelMPASParameters tlmparams {this};  // includes model parameters
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_TLMMPASPARAMETERS_H_
