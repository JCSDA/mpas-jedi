/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/Logger.h"
#include "GeometryMPAS.h"
#include "Fortran.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
GeometryMPAS::GeometryMPAS(const eckit::Configuration & conf) {
  const eckit::Configuration * configc = &conf;
  oops::Log::trace() << "============ GeometryMPAS::GeometryMPAS step 1 =============" << std::endl;
  mpas_geo_setup_f90(keyGeom_, &configc);
  oops::Log::trace() << "============ GeometryMPAS::GeometryMPAS step 2 =============" << std::endl;
}
// -----------------------------------------------------------------------------
GeometryMPAS::GeometryMPAS(const GeometryMPAS & other) {
  const int key_geo = other.keyGeom_;
  oops::Log::trace() << "============ GeometryMPAS mpas_geo_clone_f90   =============" << std::endl;
  mpas_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryMPAS::~GeometryMPAS() {
  mpas_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
void GeometryMPAS::print(std::ostream & os) const {
  int nCells;
  int nEdges;
  int nVertLevels;
  int nVertLevelsP1;
  int nCellsLocal;
  mpas_geo_info_f90(keyGeom_,nCells, nEdges, nVertLevels, nVertLevelsP1, nCellsLocal);
  os << "nCells = " << nCells << ", nEdges = " << nEdges <<", nVertLevels = "<<nVertLevels<<", nVertLevelsP1 = "<<nVertLevelsP1 << "nCellsLocal = " << nCellsLocal ;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
