/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/util/Printable.h"

#include "mpasjedi/VariableChange/Base/VariableChangeBase.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

class VariableChangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, oops::Parameters)
 public:
  // Wrapper to VariableChange parameters
  VariableChangeParametersWrapper variableChangeParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname() {return "mpas::VariableChange";}

  explicit VariableChange(const eckit::Configuration &, const Geometry &);
  ~VariableChange();

  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(State &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<VariableChangeBase> variableChange_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
