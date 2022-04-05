/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_INCREMENTMPAS_H_
#define MPASJEDI_INCREMENTMPAS_H_

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

#include "mpasjedi/IncrementMPASFortran.h"
#include "mpasjedi/IncrementMPASParameters.h"

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class GeometryMPAS;
  class StateMPAS;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class IncrementMPAS : public util::Printable,
                      public util::Serializable,
                      private util::ObjectCounter<IncrementMPAS> {
 public:
  static const std::string classname() {return "mpas::IncrementMPAS";}

  typedef IncrementMPASReadParameters ReadParameters_;
  typedef IncrementMPASWriteParameters WriteParameters_;
  typedef DiracParameters DiracParameters_;

/// Constructor, destructor
  IncrementMPAS(const GeometryMPAS &, const oops::Variables &,
                const util::DateTime &);
  IncrementMPAS(const GeometryMPAS &, const IncrementMPAS &);
  IncrementMPAS(const IncrementMPAS &, const bool);
  IncrementMPAS(const IncrementMPAS &);
  virtual ~IncrementMPAS();

/// Basic operators
  void diff(const StateMPAS &, const StateMPAS &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  IncrementMPAS & operator =(const IncrementMPAS &);
  IncrementMPAS & operator+=(const IncrementMPAS &);
  IncrementMPAS & operator-=(const IncrementMPAS &);
  IncrementMPAS & operator*=(const double &);
  void axpy(const double &, const IncrementMPAS &, const bool check = true);
  void axpy(const double &, const StateMPAS &, const bool check = true);
  double dot_product_with(const IncrementMPAS &) const;
  void schur_product_with(const IncrementMPAS &);
  void random();
  void dirac(const DiracParameters_ &);

  /// ATLAS
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);
  void getFieldSet(const oops::Variables &, atlas::FieldSet &) const;
  void getFieldSetAD(const oops::Variables &, const atlas::FieldSet &);

/// I/O and diagnostics
  void read(const ReadParameters_ &);
  void write(const WriteParameters_ &) const;
  double norm() const;

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const double &, const StateMPAS &);

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  std::shared_ptr<const GeometryMPAS> geometry() const {return geom_;}
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
  std::shared_ptr<const GeometryMPAS> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPASJEDI_INCREMENTMPAS_H_
