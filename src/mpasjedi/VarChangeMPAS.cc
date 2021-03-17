/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/Logger.h"

#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"
#include "mpasjedi/VarChangeMPAS.h"

namespace mpas {
// -----------------------------------------------------------------------------
VarChangeMPAS::VarChangeMPAS(const StateMPAS & bg,
                             const StateMPAS & fg,
                             const GeometryMPAS & resol,
                             const eckit::Configuration & config) {
  mpas_varchange_setup_f90(keyVarChange_, bg.toFortran(),
                           fg.toFortran(), resol.toFortran(),
                           config);
  oops::Log::trace() << "VarChangeMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
VarChangeMPAS::~VarChangeMPAS() {
  mpas_varchange_delete_f90(keyVarChange_);
  oops::Log::trace() << "VarChangeMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiply(const IncrementMPAS & dxa,
                             IncrementMPAS & dxm) const {
  dxm = dxa;
//  mpas_varchange_multiply_f90(keyVarChange_, dxa.toFortran(),
//                              dxm.toFortran());
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyInverse(const IncrementMPAS & dxm,
                                    IncrementMPAS & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyAD(const IncrementMPAS & dxm,
                               IncrementMPAS & dxa) const {
  dxa = dxm;
//  mpas_varchange_multiplyadjoint_f90(keyVarChange_,
//                                     dxm.toFortran(),
//                                     dxa.toFortran());
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiplyInverseAD(const IncrementMPAS & dxa,
                                      IncrementMPAS & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::print(std::ostream & os) const {
  os << "MPAS change variable";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
