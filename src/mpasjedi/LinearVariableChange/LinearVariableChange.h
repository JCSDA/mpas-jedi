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

#include "eckit/config/Configuration.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

class LinearVariableChangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, oops::Parameters)
 public:
  // Wrapper to LinearVariableChange parameters
  LinearVariableChangeParametersWrapper linearVariableChangeParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
 public:
  static const std::string classname() {return "mpas::LinearVariableChange";}

  explicit LinearVariableChange(const Geometry &, const eckit::Configuration &);
  ~LinearVariableChange();

  void changeVarTraj(const State &, const oops::Variables &);

  void changeVarTL(Increment &, const oops::Variables &) const;
  void changeVarInverseTL(Increment &, const oops::Variables &) const;
  void changeVarAD(Increment &, const oops::Variables &) const;
  void changeVarInverseAD(Increment &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  const Geometry & geom_;
  LinearVariableChangeParameters params_;
  std::unique_ptr<LinearVariableChangeBase> linearVariableChange_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
