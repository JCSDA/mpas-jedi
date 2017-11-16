/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_MPASFIELDS_H_
#define MPAS_MODEL_MPASFIELDS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "model/Geometry.h"
#include "model/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// GD add include Locations geovals
#include "/home/vagrant/jedi/code/ufo/src/ufo/GeoVaLs.h"
#include "/home/vagrant/jedi/code/ufo/src/ufo/Locations.h"


// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
}

namespace ufo {
  class Locations;
  class GeoVaLs;
}

namespace mpas {
  class Geometry;

// -----------------------------------------------------------------------------
/// Class to represent a FieldSet for the MPAS model
class Fields : public util::Printable,
                 private util::ObjectCounter<Fields> {
 public:
  static const std::string classname() {return "mpas::Fields";}

// Constructors and basic operators
  Fields(const Geometry &, const Variables &, const util::DateTime &);
  Fields(const Fields &, const Geometry &);
  Fields(const Fields &, const Variables &);
  Fields(const Fields &, const bool);
  Fields(const Fields &);
  ~Fields();

  void zero();
  void zero(const util::DateTime &);
  void dirac(const eckit::Configuration &);
  Fields & operator=(const Fields &);
  Fields & operator+=(const Fields &);
  Fields & operator-=(const Fields &);
  Fields & operator*=(const double &);
  void axpy(const double &, const Fields &);
  double dot_product_with(const Fields &) const;
  void schur_product_with(const Fields &);
  void random();

// Interpolate to given location
  void interpolate(const ufo::Locations &, ufo::GeoVaLs &) const;
  //void interpolateTL(const ufo::Locations &, ufo::GeoVaLs &) const;
  //void interpolateAD(const ufo::Locations &, const ufo::GeoVaLs &);

// Interpolate full fields
  void changeResolution(const Fields &);
  void add(const Fields &);
  void diff(const Fields &, const Fields &);

// Convert to/from generic unstructured grid
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  boost::shared_ptr<const Geometry> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

  int & toFortran() {return keyFlds_;}
  const int & toFortran() const {return keyFlds_;}

  bool isForModel(const bool) const;

 private:
  void print(std::ostream &) const;
  int keyFlds_;
  boost::shared_ptr<const Geometry> geom_;
  boost::shared_ptr<const Variables> vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_MPASFIELDS_H_
