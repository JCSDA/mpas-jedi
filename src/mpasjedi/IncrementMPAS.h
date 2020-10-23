/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_INCREMENTMPAS_H_
#define MPASJEDI_INCREMENTMPAS_H_

#include <ostream>
#include <string>
#include <memory>

#include <atlas/field.h>

#include <oops/base/GeneralizedDepartures.h>
#include <oops/base/Variables.h>
#include <oops/util/DateTime.h>
#include <oops/util/dot_product.h>
#include <oops/util/Duration.h>
#include <oops/util/ObjectCounter.h>
#include <oops/util/Printable.h>
#include <oops/util/Serializable.h>

#include "mpasjedi/IncrementMPASFortran.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class GeometryMPAS;
  class GetValuesTrajMPAS;
  class StateMPAS;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class IncrementMPAS : public oops::GeneralizedDepartures,
                       public util::Printable,
                       public util::Serializable,
                       private util::ObjectCounter<IncrementMPAS> {
 public:
  static const std::string classname() {return "mpas::IncrementMPAS";}

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
  void dirac(const eckit::Configuration &);

  /// ATLAS
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const double &, const StateMPAS &);

  std::shared_ptr<const GeometryMPAS> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  int & toFortran() {return keyInc_;}
  const int & toFortran() const {return keyInc_;}

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

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
