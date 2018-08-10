/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_VARCHANGE_H_
#define MPAS_MODEL_VARCHANGE_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class GeometryMPAS;
  class StateMPAS;
  class IncrementMPAS;

// -----------------------------------------------------------------------------
/// MPAS linear change of variable

class VarChangeMPAS: public util::Printable {
 public:
  static const std::string classname() {return "mpas::VarChangeMPAS";}

  explicit VarChangeMPAS(const eckit::Configuration &);
  ~VarChangeMPAS();

/// Set linearisation state
  void linearize(const StateMPAS &, const GeometryMPAS &);

/// Perform linear multiplications
  void multiply(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverse(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyAD(const IncrementMPAS &, IncrementMPAS &) const;
  void multiplyInverseAD(const IncrementMPAS &, IncrementMPAS &) const;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_VARCHANGE_H_
