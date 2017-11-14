/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "util/Logger.h"
#include "model/Geometry.h"
#include "model/Fortran.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & conf) {
  const eckit::Configuration * configc = &conf;
  mpas_geo_setup_f90(keyGeom_, &configc);
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) {
  const int key_geo = other.keyGeom_;
  mpas_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
  mpas_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  int nCells;
  int nEdges;
  int nVertLevels;
  int nVertLevelsP1;
  mpas_geo_info_f90(keyGeom_, nCells, nEdges, nVertLevels, nVertLevelsP1);
  os << "nCells = " << nCells << ", nVertLevels = " << nVertLevels;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
