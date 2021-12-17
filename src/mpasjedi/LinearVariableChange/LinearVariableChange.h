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

#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

class LinearVariableChangeParameters : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, oops::LinearVariableChangeParametersBase)
 public:
  // Wrapper to LinearVariableChange parameters
  LinearVariableChangeParametersWrapper linearVariableChangeParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
 public:
  static const std::string classname() {return "mpas::LinearVariableChange";}

  typedef LinearVariableChangeParameters Parameters_;

  explicit LinearVariableChange(const GeometryMPAS &, const Parameters_ &);
  ~LinearVariableChange();

  void setTrajectory(const StateMPAS &, const StateMPAS &);

  void multiply(IncrementMPAS &, const oops::Variables &) const;
  void multiplyInverse(IncrementMPAS &, const oops::Variables &) const;
  void multiplyAD(IncrementMPAS &, const oops::Variables &) const;
  void multiplyInverseAD(IncrementMPAS &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  std::shared_ptr<const GeometryMPAS> geom_;
  Parameters_ params_;
  std::unique_ptr<LinearVariableChangeBase> linearVariableChange_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
