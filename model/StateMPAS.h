/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_STATEMPAS_H_
#define MODEL_STATEMPAS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "model/StateMPASFortran.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

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
  class IncrementMPAS;

/// MPAS model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateMPAS : public util::Printable,
                  public util::Serializable,
                  private util::ObjectCounter<StateMPAS> {
 public:
  static const std::string classname() {return "mpas::StateMPAS";}

/// Constructor, destructor
  StateMPAS(const GeometryMPAS &, const oops::Variables &,
            const util::DateTime &);  // Is it used?
  StateMPAS(const GeometryMPAS &, const eckit::Configuration &);
  StateMPAS(const GeometryMPAS &, const StateMPAS &);
  StateMPAS(const StateMPAS &);
  ~StateMPAS();
//  virtual ~StateMPAS();

  StateMPAS & operator=(const StateMPAS &);
  void zero();
  void accumul(const double &, const StateMPAS &);

/// Interpolate full fields
  void changeResolution(const StateMPAS & xx);

/// Interactions with Increment
  StateMPAS & operator+=(const IncrementMPAS &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &, const GeometryMPAS &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  boost::shared_ptr<const GeometryMPAS> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  const oops::Variables & variables() const {return vars_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

 private:
  void print(std::ostream &) const override;
  F90state keyState_;
  boost::shared_ptr<const GeometryMPAS> geom_;
  oops::Variables vars_;
  util::DateTime time_;

  static const oops::Variables auxvars_;
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_STATEMPAS_H_
