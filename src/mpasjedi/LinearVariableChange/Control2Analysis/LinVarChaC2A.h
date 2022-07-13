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
#include "mpasjedi/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class Geometry;
  class State;
  class Increment;

// -----------------------------------------------------------------------------
class LinVarChaC2A : public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "mpas::LinVarChaC2A";}
  LinVarChaC2A(const State &, const State &, const Geometry &,
               const eckit::LocalConfiguration &);
  ~LinVarChaC2A();
  // Perform linear multiplications
  void changeVarTL(const Increment &, Increment &) const override;
  void changeVarInverseTL(const Increment &, Increment &) const override;
  void changeVarAD(const Increment &, Increment &) const override;
  void changeVarInverseAD(const Increment &, Increment &) const override;

 private:
  const Geometry & geom_;
  oops::Variables vars_;
  F90lvc_C2A keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
