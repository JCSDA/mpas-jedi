/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPASJEDI_INCREMENTMPASPARAMETERS_H_
#define MPASJEDI_INCREMENTMPASPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/base/WriteParametersBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

/// Parameters for Dirac
class DiracParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiracParameters, Parameters)

 public:
  // Dirac parameters
  // Number of Diracs
  oops::RequiredParameter<int> ndir{"ndir", this};
  oops::RequiredParameter<std::vector<double>> dirLats{"dirLats", this};
  oops::RequiredParameter<std::vector<double>> dirLons{"dirLons", this};
  oops::RequiredParameter<int> ildir{"ildir", this};
  oops::RequiredParameter<std::string> dirvar{"dirvar", this};
};

// -------------------------------------------------------------------------------------------------

/// Configuration options recognized by mpas_fields_mod for read_fields
class IncrementMPASReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementMPASReadParameters, Parameters)

 public:
  // Read parameters
  oops::RequiredParameter<std::string> filename{"filename", this};
  oops::RequiredParameter<util::DateTime> date{"date", this};
  oops::Parameter<std::string> stream_name{"stream name", "control", this};
};

// -------------------------------------------------------------------------------------------------

/// Configuration options recognized by mpas_fields_mod for write_fields
class IncrementMPASWriteParameters : public oops::WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(IncrementMPASWriteParameters, WriteParametersBase)

 public:
  // Write parameters
  oops::RequiredParameter<std::string> filename{"filename", this};
  oops::Parameter<std::string> stream_name{"stream name", "da_state", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPASJEDI_INCREMENTMPASPARAMETERS_H_
