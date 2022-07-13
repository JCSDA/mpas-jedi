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

  void changeVarTL(const Increment &, Increment &) const override;
  void changeVarInverseTL(const Increment &, Increment &) const override;
  void changeVarAD(const Increment &, Increment &) const override;
  void changeVarInverseAD(const Increment &, Increment &) const override;

 private:
  const Geometry & geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
