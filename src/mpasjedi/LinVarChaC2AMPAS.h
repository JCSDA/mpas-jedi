/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_LINVARCHAC2AMPAS_H_
#define MODEL_LINVARCHAC2AMPAS_H_

#include <ostream>
#include <string>

#include <eckit/config/Configuration.h>

#include <oops/base/Variables.h>
#include <oops/util/Printable.h>

#include "mpasjedi/Fortran.h"
#include "mpasjedi/GeometryMPAS.h"

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
/// MPAS linear change of variable

class LinVarChaC2AMPAS: public util::Printable {
 public:
  static const std::string classname() {return "mpas::LinVarChaC2AMPAS";}
  explicit LinVarChaC2AMPAS(const StateMPAS &, const StateMPAS &,
                         const GeometryMPAS &, const eckit::Configuration &);
  ~LinVarChaC2AMPAS();
  // Perform linear multiplications
  void multiply(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverse(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyAD(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverseAD(const IncrementMPAS &, IncrementMPAS &) const;

 private:
  oops::Variables vars_;
  F90lvcc2a keyLinVarChaC2A_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_LINVARCHAC2AMPAS_H_
