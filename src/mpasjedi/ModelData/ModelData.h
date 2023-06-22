/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>
#include "eckit/config/LocalConfiguration.h"
#include "oops/util/Printable.h"

namespace mpas {
  class Geometry;
}

// -------------------------------------------------------------------------------------------------

namespace mpas {

class ModelData : public util::Printable {
 public:
  static const std::string classname() {return "mpas::ModelData";}

  explicit ModelData(const Geometry &) {}
  ~ModelData() {}

  const eckit::LocalConfiguration modelData() const {return eckit::LocalConfiguration();}

 private:
  void print(std::ostream & os) const {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
