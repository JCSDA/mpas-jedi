/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_MODELMPAS_H_
#define MODEL_MODELMPAS_H_

#include <ostream>
#include <string>

#include "model/GeometryMPAS.h"
#include "model/MPASTraits.h"

#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace mpas {
  class ModelBiasMPAS;
  class StateMPAS;

// -----------------------------------------------------------------------------
/// MPAS model definition.
/*!
 *  MPAS nonlinear model definition and configuration parameters.
 */

class ModelMPAS: public oops::ModelBase<MPASTraits>,
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
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyModel_;
  util::Duration tstep_;
  const GeometryMPAS geom_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_MODELMPAS_H_
