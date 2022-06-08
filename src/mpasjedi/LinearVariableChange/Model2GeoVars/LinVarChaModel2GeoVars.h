/*
 * (C) Copyright 2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "mpasjedi/LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.interface.h"
#include "mpasjedi/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class Geometry;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVars : public LinearVariableChangeBase,
                               private util::ObjectCounter<LinVarChaModel2GeoVars> {
 public:
  static const std::string classname() {return "mpas::LinVarChaModel2GeoVars";}

  explicit LinVarChaModel2GeoVars(const State &, const State &, const Geometry &,
                         const eckit::LocalConfiguration &);
  ~LinVarChaModel2GeoVars();

  void multiply(const Increment &, Increment &) const override;
  void multiplyInverse(const Increment &, Increment &) const override;
  void multiplyAD(const Increment &, Increment &) const override;
  void multiplyInverseAD(const Increment &, Increment &) const override;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
