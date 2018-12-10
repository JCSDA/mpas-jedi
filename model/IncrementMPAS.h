/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_INCREMENTMPAS_H_
#define MODEL_INCREMENTMPAS_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "FieldsMPAS.h"
#include "GeometryMPAS.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/dot_product.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
  class UnstructuredGrid;
}

namespace mpas {
  class ModelBiasIncrementMPAS;
  class ErrorCovarianceMPAS;
  class StateMPAS;
  class GetValuesTrajMPAS;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class IncrementMPAS : public oops::GeneralizedDepartures,
                       public util::Printable,
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
  IncrementMPAS & operator =(const IncrementMPAS &);
  IncrementMPAS & operator+=(const IncrementMPAS &);
  IncrementMPAS & operator-=(const IncrementMPAS &);
  IncrementMPAS & operator*=(const double &);
  void axpy(const double &, const IncrementMPAS &, const bool check = true);
  double dot_product_with(const IncrementMPAS &) const;
  void schur_product_with(const IncrementMPAS &);
  void random();
  void dirac(const eckit::Configuration &);

/// Get increment values at observation locations
  void getValuesTL(const ufo::Locations &, const oops::Variables &,
                   ufo::GeoVaLs &, const GetValuesTrajMPAS &) const;
  void getValuesAD(const ufo::Locations &, const oops::Variables &,
                   const ufo::GeoVaLs &, const GetValuesTrajMPAS &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// Unstructured grid
  void ug_coord(oops::UnstructuredGrid &, const int &) const;
  void field_to_ug(oops::UnstructuredGrid &, const int &) const;
  void field_from_ug(const oops::UnstructuredGrid &);

/// Access to fields
  FieldsMPAS & fields() {return *fields_;}
  const FieldsMPAS & fields() const {return *fields_;}

  boost::shared_ptr<const GeometryMPAS> geometry() const {
    return fields_->geometry();
  }

/// Other
  void accumul(const double &, const StateMPAS &);

/// Data
 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<FieldsMPAS> fields_;
  boost::scoped_ptr<FieldsMPAS> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_INCREMENTMPAS_H_
