/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "mpasjedi/LinearVariableChange/Control2Analysis/LinVarChaC2A.interface.h"
#include "mpasjedi/MPASTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class GeometryMPAS;
  class StateMPAS;
  class IncrementMPAS;

// -----------------------------------------------------------------------------
class LinVarChaC2A : public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "mpas::LinVarChaC2A";}
  LinVarChaC2A(const StateMPAS &, const StateMPAS &, const GeometryMPAS &,
               const eckit::LocalConfiguration &);
  ~LinVarChaC2A();
  // Perform linear multiplications
  void multiply(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyInverse(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyAD(const IncrementMPAS &, IncrementMPAS &) const override;
  void multiplyInverseAD(const IncrementMPAS &, IncrementMPAS &) const override;

 private:
  std::shared_ptr<const GeometryMPAS> geom_;
  oops::Variables vars_;
  F90lvc_C2A keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
