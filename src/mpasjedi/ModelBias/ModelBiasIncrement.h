/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <iostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace mpas {
  class ModelBias;
  class ModelBiasCovariance;
  class Geometry;


// -----------------------------------------------------------------------------

class ModelBiasIncrement : public util::Printable,
                               public util::Serializable {
 public:
/// Constructor, destructor
  ModelBiasIncrement(const Geometry &, const eckit::Configuration &) {}
  ModelBiasIncrement(const ModelBiasIncrement &, const bool) {}
  ModelBiasIncrement(const ModelBiasIncrement &,
                         const eckit::Configuration &) {}
  ~ModelBiasIncrement() {}

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &) {}
  void zero() {}
  ModelBiasIncrement & operator=(const ModelBiasIncrement &)
                                    {return *this;}
  ModelBiasIncrement & operator+=(const ModelBiasIncrement &)
                                     {return *this;}
  ModelBiasIncrement & operator-=(const ModelBiasIncrement &)
                                     {return *this;}
  ModelBiasIncrement & operator*=(const real_type) {return *this;}
  void axpy(const real_type, const ModelBiasIncrement &) {}
  real_type dot_product_with(const ModelBiasIncrement &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  real_type norm() const {return 0.0;}

/// Serialization
  size_t serialSize() const override { return 0; }
  void serialize(std::vector<real_type> &) const override {}
  void deserialize(const std::vector<real_type> &, size_t &) override {}

 private:
  explicit ModelBiasIncrement(const ModelBiasCovariance &);
  void print(std::ostream & os) const override {}
};

// -----------------------------------------------------------------------------

}  // namespace mpas
