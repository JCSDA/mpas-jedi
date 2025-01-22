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

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "mpasjedi/State/State.interface.h"
#include "mpasjedi/State/StateParameters.h"

namespace ufo {
  class GeoVaLs;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class Geometry;
  class Increment;

/// MPAS model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class State : public util::Printable,
                  public util::Serializable,
                  private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "mpas::State";}

/// Constructor, destructor
  State(const Geometry &, const oops::Variables &,
            const util::DateTime &);
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const oops::Variables &, const State &);
  State(const State &);
  ~State();
//  virtual ~State();

  State & operator=(const State &);
  void zero();
  void accumul(const real_type &, const State &);

/// Interpolate full fields
  void changeResolution(const State & xx);

/// Interactions with Increment
  State & operator+=(const Increment &);

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<real_type> &) const override;
  void deserialize(const std::vector<real_type> &, size_t &) override;
  void transpose(const State & DistState, const eckit::mpi::Comm & global,
     const int ensNum, const int transNum) {
     throw eckit::NotImplemented("MPASJEDI State::transpose not implemented", Here());
  }

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  real_type norm() const;

  const Geometry & geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  const oops::Variables & variables() const {return vars_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  // Accessors to the ATLAS fieldset
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

 private:
/// Analytic intialization
  void analytic_init(const eckit::Configuration &);

  void print(std::ostream &) const override;
  F90state keyState_;
  const Geometry & geom_;
  oops::Variables vars_;
  util::DateTime time_;

  oops::Variables stateVars();
};

// -----------------------------------------------------------------------------

}  // namespace mpas
