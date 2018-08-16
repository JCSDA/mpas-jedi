/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_FIELDSMPAS_H_
#define MPAS_MODEL_FIELDSMPAS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "GeometryMPAS.h"
#include "GetValuesTrajMPAS.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
//  class Variables;
}

namespace ufo {
  class GeoVaLs;
}

namespace ioda {
  class Locations;
}

namespace mpas {
// -----------------------------------------------------------------------------
/// Class to represent a FieldSet for the MPAS model
class FieldsMPAS : public util::Printable,
                    private util::ObjectCounter<FieldsMPAS> {
 public:
  static const std::string classname() {return "mpas::FieldsMPAS";}

// Constructors and basic operators
  FieldsMPAS(const GeometryMPAS &, const oops::Variables &,
             const util::DateTime &);
  FieldsMPAS(const FieldsMPAS &, const GeometryMPAS &);
  FieldsMPAS(const FieldsMPAS &, const oops::Variables &);
  FieldsMPAS(const FieldsMPAS &, const bool);
  FieldsMPAS(const FieldsMPAS &);
  ~FieldsMPAS();

  void zero();
  void zero(const util::DateTime &);
  FieldsMPAS & operator=(const FieldsMPAS &);
  FieldsMPAS & operator+=(const FieldsMPAS &);
  FieldsMPAS & operator-=(const FieldsMPAS &);
  FieldsMPAS & operator*=(const double &);
  void axpy(const double &, const FieldsMPAS &);
  double dot_product_with(const FieldsMPAS &) const;
  void schur_product_with(const FieldsMPAS &);
  void random();
  void dirac(const eckit::Configuration &);

/// Get state values or increments at observation locations
  void getValues(const ioda::Locations &, const oops::Variables &,
                 ufo::GeoVaLs &) const;
  void getValues(const ioda::Locations &, const oops::Variables &,
                 ufo::GeoVaLs &, const GetValuesTrajMPAS &) const;
  void getValuesTL(const ioda::Locations &, const oops::Variables &,
                   ufo::GeoVaLs &, const GetValuesTrajMPAS &) const;
  void getValuesAD(const ioda::Locations &, const oops::Variables &,
                   const ufo::GeoVaLs &, const GetValuesTrajMPAS &);

// Interpolate full fields
  void changeResolution(const FieldsMPAS &);
  void add(const FieldsMPAS &);
  void diff(const FieldsMPAS &, const FieldsMPAS &);

// Unstructured grid
  void ug_coord(oops::UnstructuredGrid &, const int &) const;
  void field_to_ug(oops::UnstructuredGrid &, const int &) const;
  void field_from_ug(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &, const GeometryMPAS &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  boost::shared_ptr<const GeometryMPAS> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

  int & toFortran() {return keyFlds_;}
  const int & toFortran() const {return keyFlds_;}

 private:
  void print(std::ostream &) const;
  F90flds keyFlds_;
  boost::shared_ptr<const GeometryMPAS> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_FIELDSMPAS_H_
