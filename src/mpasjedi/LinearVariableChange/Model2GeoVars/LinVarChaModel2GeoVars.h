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
#include "mpasjedi/MPASTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class GeometryMPAS;
  class IncrementMPAS;
  class StateMPAS;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVars : public LinearVariableChangeBase,
                               private util::ObjectCounter<LinVarChaModel2GeoVars> {
 public:
  static const std::string classname() {return "mpas::LinVarChaModel2GeoVars";}

  explicit LinVarChaModel2GeoVars(const StateMPAS &, const StateMPAS &, const GeometryMPAS &,
                         const eckit::LocalConfiguration &);
  ~LinVarChaModel2GeoVars();

  void multiply(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyInverse(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyAD(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyInverseAD(const IncrementMPAS &, IncrementMPAS &) const override;

 private:
  std::shared_ptr<const GeometryMPAS> geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
