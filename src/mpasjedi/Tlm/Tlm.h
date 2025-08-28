/*
 * (C) Copyright 2017-2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "mpasjedi/Fortran.h"
#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Model/Model.h"

// Forward declarations

namespace mpas {
  class State;
  class Increment;
  class ModelBias;
  class ModelBiasIncrement;

// -----------------------------------------------------------------------------
/// MPAS linear model definition.
/*!
 *  MPAS linear model definition and configuration parameters.
 */

class Tlm: public util::Printable,
           private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "mpas::Tlm";}
  static std::vector<std::string> names() {return {"MPASTLM"};}

  Tlm(const Geometry &, const eckit::Configuration &);
  ~Tlm();

/// Model trajectory computation
  void setTrajectory(const State &, State &, const ModelBias &);

/// Run TLM and its adjoint
  void initializeTL(Increment &) const;
  void stepTL(Increment &, const ModelBiasIncrement &) const;
  void finalizeTL(Increment &) const;

  void initializeAD(Increment &) const;
  void stepAD(Increment &, ModelBiasIncrement &) const;
  void finalizeAD(Increment &) const;

/// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return tstep_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry resol_;
  std::map< util::DateTime, F90traj> traj_;
  const Model lrmodel_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
