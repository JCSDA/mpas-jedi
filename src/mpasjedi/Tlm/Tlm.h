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

#include "oops/base/Variables.h"
#include "oops/interface/LinearModelBase.h"

#include "mpasjedi/Fortran.h"
#include "mpasjedi/Model/Model.h"
#include "mpasjedi/Tlm/TlmParameters.h"
#include "mpasjedi/Traits.h"

// Forward declarations

namespace mpas {
  class State;
  class Increment;
  // class ModelBias;
  // class ModelBiasIncrement;

// -----------------------------------------------------------------------------
/// MPAS linear model definition.
/*!
 *  MPAS linear model definition and configuration parameters.
 */

class Tlm: public oops::interface::LinearModelBase<Traits>,
                private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "mpas::Tlm";}

  Tlm(const Geometry &, const eckit::Configuration &);
  ~Tlm();

/// Model trajectory computation
  void setTrajectory(const State &, State &, const ModelBias &)
                    override;

/// Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &) const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &) const override;
  void finalizeAD(Increment &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Geometry & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

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
  oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
