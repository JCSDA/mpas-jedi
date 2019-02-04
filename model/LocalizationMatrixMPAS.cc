/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/LocalizationMatrixMPAS.h"

#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
LocalizationMatrixMPAS::LocalizationMatrixMPAS(const GeometryMPAS & resol,
                                          const eckit::Configuration & config) {
//  const eckit::Configuration * configc = &config;
}
// -----------------------------------------------------------------------------
LocalizationMatrixMPAS::~LocalizationMatrixMPAS() {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixMPAS::multiply(IncrementMPAS & dx) const {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixMPAS::print(std::ostream & os) const {
  os << "LocalizationMatrixMPAS::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace mpas
