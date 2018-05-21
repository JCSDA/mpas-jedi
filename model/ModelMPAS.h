/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_MODELMPAS_H_
#define MPAS_MODEL_MODELMPAS_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "Fortran.h"
#include "GeometryMPAS.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class ModelBiasMPAS;
  class FieldsMPAS;
  class StateMPAS;

// -----------------------------------------------------------------------------
/// MPAS model definition.
/*!
 *  MPAS nonlinear model definition and configuration parameters.
 */

class ModelMPAS: public util::Printable,
               private boost::noncopyable,
               private util::ObjectCounter<ModelMPAS> {
 public:
  static const std::string classname() {return "mpas::ModelMPAS";}

  ModelMPAS(const GeometryMPAS &, const eckit::Configuration &);
  ~ModelMPAS();

/// Prepare model integration
  void initialize(StateMPAS &) const;

/// Model integration
  void step(StateMPAS &, const ModelBiasMPAS &) const;
  int saveTrajectory(StateMPAS &, const ModelBiasMPAS &) const;

/// Finish model integration
  void finalize(StateMPAS &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryMPAS geom_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_MODELMPAS_H_
