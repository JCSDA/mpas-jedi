/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_STATEMPAS_H_
#define MPAS_MODEL_STATEMPAS_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "FieldsMPAS.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

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
  class IncrementMPAS;

/// MPAS model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateMPAS : public util::Printable,
                private util::ObjectCounter<StateMPAS> {
 public:
  static const std::string classname() {return "mpas::StateMPAS";}

/// Constructor, destructor
  StateMPAS(const GeometryMPAS &, const oops::Variables &, const util::DateTime &);  // Is it used?
  StateMPAS(const GeometryMPAS &, const eckit::Configuration &);
  StateMPAS(const GeometryMPAS &, const StateMPAS &);
  StateMPAS(const StateMPAS &);
  virtual ~StateMPAS();
  StateMPAS & operator=(const StateMPAS &);

/// Interpolate to observation location
  void interpolate(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;

/// Interpolate full fields
  void changeResolution(const StateMPAS & xx);

/// Interactions with Increment
  StateMPAS & operator+=(const IncrementMPAS &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}

/// Access to fields
  FieldsMPAS & fields() {return *fields_;}
  const FieldsMPAS & fields() const {return *fields_;}

  boost::shared_ptr<const GeometryMPAS> geometry() const {
    return fields_->geometry();
  }

/// Other
  void activateModel();
  void deactivateModel();

  void zero();
  void accumul(const double &, const StateMPAS &);

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<FieldsMPAS> fields_;
  boost::scoped_ptr<FieldsMPAS> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPAS_MODEL_STATEMPAS_H_
