/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
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
  oops::RequiredParameter<std::vector<real_type>> dirLats{"dirLats", this};
  oops::RequiredParameter<std::vector<real_type>> dirLons{"dirLons", this};
  oops::RequiredParameter<int> ildir{"ildir", this};
  oops::RequiredParameter<std::string> dirvar{"dirvar", this};
};

// -------------------------------------------------------------------------------------------------

/// Configuration options recognized by mpas_fields_mod for read_fields
class IncrementReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementReadParameters, Parameters)

 public:
  // Read parameters
  oops::RequiredParameter<std::string> filename{"filename", this};
  oops::OptionalParameter<util::DateTime> date{"date", this};
  oops::Parameter<std::string> stream_name{"stream name", "control", this};
};

// -------------------------------------------------------------------------------------------------

/// Configuration options recognized by mpas_fields_mod for write_fields
class IncrementWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementWriteParameters, Parameters)

 public:
  // Write parameters
  oops::OptionalParameter<std::string> type{"type", this};
  oops::OptionalParameter<std::string> exp{"exp", this};
  oops::OptionalParameter<int> member{"member", this};
  oops::OptionalParameter<std::string> memberPattern{"member pattern", this};
  oops::OptionalParameter<util::DateTime> date{"date", this};
  oops::OptionalParameter<int> iteration{"iteration", this};
  oops::OptionalParameter<std::string> prefix{"prefix", this};
  oops::Parameter<bool> dateCols{"date colons", true, this};
  oops::OptionalParameter<std::string> filename{"filename", this};
  oops::Parameter<std::string> stream_name{"stream name", "da_state", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
