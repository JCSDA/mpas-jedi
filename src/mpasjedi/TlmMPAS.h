/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_TLMMPAS_H_
#define MPASJEDI_TLMMPAS_H_

#include <map>
#include <ostream>
#include <string>

#include <oops/base/LinearModelBase.h>
#include <oops/base/Variables.h>

#include "mpasjedi/Fortran.h"
#include "mpasjedi/ModelMPAS.h"
#include "mpasjedi/MPASTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class StateMPAS;
  class IncrementMPAS;

// -----------------------------------------------------------------------------
/// LFRic linear model definition.
/*!
 *  LFRic linear model definition and configuration parameters.
 */

class TlmMPAS: public oops::LinearModelBase<MPASTraits>,
                private util::ObjectCounter<TlmMPAS> {
 public:
  static const std::string classname() {return "mpas::TlmMPAS";}

  TlmMPAS(const GeometryMPAS &, const eckit::Configuration &);
  ~TlmMPAS();

/// Model trajectory computation
  void setTrajectory(const StateMPAS &, StateMPAS &, const ModelBiasMPAS &)
                    override;

/// Run TLM and its adjoint
  void initializeTL(IncrementMPAS &) const override;
  void stepTL(IncrementMPAS &, const ModelBiasIncrementMPAS &) const override;
  void finalizeTL(IncrementMPAS &) const override;

  void initializeAD(IncrementMPAS &) const override;
  void stepAD(IncrementMPAS &, ModelBiasIncrementMPAS &) const override;
  void finalizeAD(IncrementMPAS &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const GeometryMPAS & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryMPAS resol_;
  std::map< util::DateTime, F90traj> traj_;
  const ModelMPAS lrmodel_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_TLMMPAS_H_
