/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_TLMIDMPAS_H_
#define MODEL_TLMIDMPAS_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/LinearModelBase.h"
#include "oops/base/Variables.h"

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "MPASTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {
  class StateMPAS;
  class IncrementMPAS;

// -----------------------------------------------------------------------------
/// MPAS linear identity model definition.
/*!
 *  MPAS linear identity model definition and configuration parameters.
 */

class TlmIdMPAS: public oops::LinearModelBase<MPASTraits>,
                  private util::ObjectCounter<TlmIdMPAS> {
 public:
  static const std::string classname() {return "mpas::TlmIdMPAS";}

  TlmIdMPAS(const GeometryMPAS &, const eckit::Configuration &);
  ~TlmIdMPAS();

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

// Data
  int keyConfig_;
  util::Duration tstep_;
  const GeometryMPAS resol_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_TLMIDMPAS_H_
