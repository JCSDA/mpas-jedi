/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_MODELBIASMPAS_H_
#define MODEL_MODELBIASMPAS_H_

#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mpas {
  class GeometryMPAS;
  class ModelBiasIncrementMPAS;

/// Model error for the MPAS model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 * This is sometimes referred to as augmented state or augmented
 * control variable in the litterature.
 * The augmented state is understood here as an augmented 4D state.
 */

// -----------------------------------------------------------------------------

class ModelBiasMPAS : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ModelBiasMPAS> {
 public:
  static const std::string classname() {return "mpas::ModelBiasMPAS";}

  ModelBiasMPAS(const GeometryMPAS &, const eckit::Configuration &) {}
  ModelBiasMPAS(const GeometryMPAS &, const ModelBiasMPAS &) {}
  ModelBiasMPAS(const ModelBiasMPAS &, const bool) {}
  ~ModelBiasMPAS() {}

  ModelBiasMPAS & operator+=(const ModelBiasIncrementMPAS &) {return *this;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_MODELBIASMPAS_H_
