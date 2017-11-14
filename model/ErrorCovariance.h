/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef MPAS_MODEL_MPASERRORCOVARIANCE_H_
#define MPAS_MODEL_MPASERRORCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "model/Fortran.h"
#include "model/Geometry.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace mpas {
  class Increment;
  class State;
  class Variables;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for MPAS model.

class ErrorCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ErrorCovariance> {
 public:
  static const std::string classname() {return "mpas::ErrorCovariance";}

  ErrorCovariance(const Geometry &, const Variables &, const eckit::Configuration &, const State &) {}
  ~ErrorCovariance() {}

  void linearize(const State &, const Geometry &) {}
  void multiply(const Increment &, Increment &) const {}
  void inverseMultiply(const Increment &, Increment &) const {}
  void randomize(Increment &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_MPASERRORCOVARIANCE_H_
