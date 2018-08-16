/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_ERRORCOVARIANCEMPAS_H_
#define MPAS_MODEL_ERRORCOVARIANCEMPAS_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "Fortran.h"
#include "GeometryMPAS.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace mpas {
  class IncrementMPAS;
  class StateMPAS;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for LFRic

class ErrorCovarianceMPAS : public util::Printable,
                             private boost::noncopyable,
                             private util::ObjectCounter<ErrorCovarianceMPAS> {
 public:
  static const std::string classname() {return "mpas::ErrorCovarianceMPAS";}

  ErrorCovarianceMPAS(const GeometryMPAS &, const oops::Variables &,
                      const eckit::Configuration &, const StateMPAS &,
                      const StateMPAS &);
  ~ErrorCovarianceMPAS();

  void linearize(const StateMPAS &, const GeometryMPAS &);
  void multiply(const IncrementMPAS &, IncrementMPAS &) const;
  void inverseMultiply(const IncrementMPAS &, IncrementMPAS &) const;
  void randomize(IncrementMPAS &) const;

 private:
  void print(std::ostream &) const;
  F90bmat keyFtnConfig_;
  boost::scoped_ptr<const GeometryMPAS> geom_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_ERRORCOVARIANCEMPAS_H_
