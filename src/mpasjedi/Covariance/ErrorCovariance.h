/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/Covariance/ErrorCovariance.interface.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace mpas {
  class Geometry;
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for LFRic

class ErrorCovariance : public util::Printable,
                             private eckit::NonCopyable,
                             private util::ObjectCounter<ErrorCovariance> {
 public:
  static const std::string classname() {return "mpas::ErrorCovariance";}

  ErrorCovariance(const Geometry &, const oops::Variables &,
                      const eckit::Configuration &, const State &,
                      const State &);
  ~ErrorCovariance();

  void linearize(const State &, const Geometry &);
  void multiply(const Increment &, Increment &) const;
  void inverseMultiply(const Increment &, Increment &) const;
  void randomize(Increment &) const;

 private:
  void print(std::ostream &) const;
  F90bmat keyErrCov_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
