/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_MODELBIASINCREMENTMPAS_H_
#define MODEL_MODELBIASINCREMENTMPAS_H_

#include <iostream>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mpas {
  class ModelBiasMPAS;
  class ModelBiasCovarianceMPAS;
  class GeometryMPAS;

// -----------------------------------------------------------------------------

class ModelBiasIncrementMPAS : public util::Printable {
 public:
/// Constructor, destructor
  ModelBiasIncrementMPAS(const GeometryMPAS &, const eckit::Configuration &) {}
  ModelBiasIncrementMPAS(const ModelBiasIncrementMPAS &, const bool) {}
  ModelBiasIncrementMPAS(const ModelBiasIncrementMPAS &,
                         const eckit::Configuration &) {}
  ~ModelBiasIncrementMPAS() {}

/// Linear algebra operators
  void diff(const ModelBiasMPAS &, const ModelBiasMPAS &) {}
  void zero() {}
  ModelBiasIncrementMPAS & operator=(const ModelBiasIncrementMPAS &)
                                    {return *this;}
  ModelBiasIncrementMPAS & operator+=(const ModelBiasIncrementMPAS &)
                                     {return *this;}
  ModelBiasIncrementMPAS & operator-=(const ModelBiasIncrementMPAS &)
                                     {return *this;}
  ModelBiasIncrementMPAS & operator*=(const double) {return *this;}
  void axpy(const double, const ModelBiasIncrementMPAS &) {}
  double dot_product_with(const ModelBiasIncrementMPAS &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  explicit ModelBiasIncrementMPAS(const ModelBiasCovarianceMPAS &);
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_MODELBIASINCREMENTMPAS_H_
