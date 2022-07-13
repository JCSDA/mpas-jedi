/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "mpasjedi/Traits.h"
#include "mpasjedi/VariableChange/Base/VariableChangeBase.h"
#include "mpasjedi/VariableChange/Model2GeoVars/VarChaModel2GeoVars.interface.h"

namespace mpas {
  class State;
  class Geometry;

// -------------------------------------------------------------------------------------------------

class VarChaModel2GeoVars : public VariableChangeBase,
                            private util::ObjectCounter<VarChaModel2GeoVars> {
 public:
  static const std::string classname() {return "mpas::VarChaModel2GeoVars";}
  VarChaModel2GeoVars(const Geometry &, const eckit::LocalConfiguration &);
  ~VarChaModel2GeoVars();
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  F90vc_M2G keyFtnConfig_;
  const Geometry & geom_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
