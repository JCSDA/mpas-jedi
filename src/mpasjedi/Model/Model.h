/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Model/Model.interface.h"

// Forward declarations

namespace oops {
  class Variables;
}

namespace mpas {
  class ModelBias;
  class ModelParameters;
  class State;

// -----------------------------------------------------------------------------
/// MPAS model definition.
/*!
 *  MPAS nonlinear model definition and configuration parameters.
 */

class Model: public util::Printable,
             private util::ObjectCounter<Model> {
 public:
  static const std::string classname() {return "mpas::Model";}
  static std::vector<std::string> names() {return {"MPAS"};}

  Model(const Geometry &, const eckit::Configuration &);
  ~Model();

/// Prepare model integration
  void initialize(State &) const;

/// Model integration
  void step(State &, const ModelBias &) const;
  int saveTrajectory(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const override;
  F90model keyModel_;
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
