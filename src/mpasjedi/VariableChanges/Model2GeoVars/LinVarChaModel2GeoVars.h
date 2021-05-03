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

#include "oops/util/Printable.h"

#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/VariableChanges/Model2GeoVars/LinVarChaModel2GeoVars.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class GeometryMPAS;
  class StateMPAS;
  class IncrementMPAS;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVars : public util::Printable,
                               private util::ObjectCounter<LinVarChaModel2GeoVars> {
 public:
  static const std::string classname() {return "mpas::LinVarChaModel2GeoVars";}

  explicit LinVarChaModel2GeoVars(const StateMPAS &, const StateMPAS &, const GeometryMPAS &,
                                  const eckit::Configuration &);
  ~LinVarChaModel2GeoVars();

  void multiply(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverse(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyAD(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverseAD(const IncrementMPAS &, IncrementMPAS &) const;

 private:
  std::shared_ptr<const GeometryMPAS> geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------

}  // namespace mpas
