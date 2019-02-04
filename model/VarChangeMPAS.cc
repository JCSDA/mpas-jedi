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
#include "oops/util/Logger.h"

namespace mpas {
// -----------------------------------------------------------------------------
VarChangeMPAS::VarChangeMPAS(const StateMPAS & bg,
                             const StateMPAS & fg,
                             const GeometryMPAS & resol,
                             const eckit::Configuration & conf) {
    const eckit::Configuration * configc = &conf;
    mpas_varchange_setup_f90(keyFtnConfig_, bg.toFortran(),
                             fg.toFortran(), resol.toFortran(),
                             &configc);
    oops::Log::trace() << "VarChangeMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
VarChangeMPAS::~VarChangeMPAS() {
    mpas_varchange_delete_f90(keyFtnConfig_);
    oops::Log::trace() << "VarChangeMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChangeMPAS::multiply(const IncrementMPAS & dxa,
                             IncrementMPAS & dxm) const {
  dxm = dxa;
//  mpas_varchange_multiply_f90(keyFtnConfig_, dxa.toFortran(),
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
//  mpas_varchange_multiplyadjoint_f90(keyFtnConfig_,
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
