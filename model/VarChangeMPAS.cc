/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/VarChangeMPAS.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryMPAS.h"
#include "model/IncrementMPAS.h"
#include "model/StateMPAS.h"

namespace mpas {
// -----------------------------------------------------------------------------
VarChangeMPAS::VarChangeMPAS(const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
VarChangeMPAS::~VarChangeMPAS() {}
// -----------------------------------------------------------------------------
void VarChangeMPAS::linearize(const StateMPAS &, const GeometryMPAS &) {}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiply(const IncrementMPAS & dxa, IncrementMPAS & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyInverse(const IncrementMPAS & dxm, IncrementMPAS & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyAD(const IncrementMPAS & dxm, IncrementMPAS & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyInverseAD(const IncrementMPAS & dxa, IncrementMPAS & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::print(std::ostream & os) const {
  os << "MPAS change variable";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
