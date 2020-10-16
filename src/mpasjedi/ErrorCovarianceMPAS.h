/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_ERRORCOVARIANCEMPAS_H_
#define MPASJEDI_ERRORCOVARIANCEMPAS_H_

#include <memory>
#include <ostream>
#include <string>

#include <eckit/config/Configuration.h>
#include <eckit/memory/NonCopyable.h>

#include <oops/util/DateTime.h>
#include <oops/util/ObjectCounter.h>
#include <oops/util/Printable.h>

#include "mpasjedi/Fortran.h"
#include "mpasjedi/GeometryMPAS.h"

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
                             private eckit::NonCopyable,
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
  F90bmat keyErrCov_;
  std::unique_ptr<const GeometryMPAS> geom_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_ERRORCOVARIANCEMPAS_H_
