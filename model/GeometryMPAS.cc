/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "util/Logger.h"
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
  mpas_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryMPAS::~GeometryMPAS() {
  mpas_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
void GeometryMPAS::print(std::ostream & os) const {
  int nx;
  int ny;
  mpas_geo_info_f90(keyGeom_, nx, ny);
  os << "nx = " << nx << ", ny = " << ny;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
