/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/LocalIncrement.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "mpasjedi/Increment/Increment.interface.h"
#include "mpasjedi/Increment/IncrementParameters.h"

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class Geometry;
  class State;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class Increment : public util::Printable,
                      public util::Serializable,
                      private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "mpas::Increment";}

  typedef IncrementReadParameters ReadParameters_;
  typedef IncrementWriteParameters WriteParameters_;
  typedef DiracParameters DiracParameters_;

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &,
                const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  Increment(const Increment &);
  virtual ~Increment();

/// Basic operators
  void diff(const State &, const State &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const real_type &);
  void axpy(const real_type &, const Increment &, const bool check = true);
  void axpy(const real_type &, const State &, const bool check = true);
  real_type dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void dirac(const DiracParameters_ &);
  std::vector<double> rmsByLevel(const std::string &) const;

  void toFieldSet(atlas::FieldSet &) const;
  void toFieldSetAD(const atlas::FieldSet &);
  void fromFieldSet(const atlas::FieldSet &);

/// I/O and diagnostics
  void read(const ReadParameters_ &);
  void write(const WriteParameters_ &) const;
  real_type norm() const;

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const real_type &, const State &);

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<real_type> &) const override;
  void deserialize(const std::vector<real_type> &, size_t &) override;

  const Geometry & geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  int & toFortran() {return keyInc_;}
  const int & toFortran() const {return keyInc_;}

/// Data
 private:
  void print(std::ostream &) const override;
  F90inc keyInc_;
  const Geometry & geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
