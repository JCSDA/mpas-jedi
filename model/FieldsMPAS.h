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
//#include "oops/base/Variables.h"
#include "VariablesMPAS.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
  class Variables;
}

namespace ufo {
  class Locations;
  class GeoVaLs;
}

namespace mpas {
// -----------------------------------------------------------------------------
/// Class to represent a FieldSet for the MPAS model
class FieldsMPAS : public util::Printable,
                    private util::ObjectCounter<FieldsMPAS> {
 public:
  static const std::string classname() {return "mpas::FieldsMPAS";}

// Constructors and basic operators
  FieldsMPAS(const GeometryMPAS &, const oops::Variables &, const util::DateTime &);
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

// Interpolate to given location
  void interpolate(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
  void interpolateTL(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
  void interpolateAD(const ufo::Locations &, const oops::Variables &, const ufo::GeoVaLs &);

// Interpolate full fields
  void changeResolution(const FieldsMPAS &);
  void add(const FieldsMPAS &);
  void diff(const FieldsMPAS &, const FieldsMPAS &);

// Define and convert to/from unstructured grid
  void define(oops::UnstructuredGrid &) const;
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
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
//  oops::Variables vars_;
  const VariablesMPAS vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_FIELDSMPAS_H_
