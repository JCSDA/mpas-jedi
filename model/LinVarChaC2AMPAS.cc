/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/LinVarChaC2AMPAS.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryMPAS.h"
#include "model/IncrementMPAS.h"
#include "model/StateMPAS.h"
#include "oops/util/Logger.h"
#include "Fortran.h"

namespace mpas {
// -----------------------------------------------------------------------------
LinVarChaC2AMPAS::LinVarChaC2AMPAS(const StateMPAS & bg,
                             const StateMPAS & fg,
                             const GeometryMPAS & resol,
                             const eckit::Configuration & config) {
  mpas_linvarcha_c2a_setup_f90(keyLinVarChaC2A_, bg.toFortran(),
                               fg.toFortran(), resol.toFortran(),
                               config);
  oops::Log::trace() << "LinVarChaC2AMPAS created" << std::endl;
}
// -----------------------------------------------------------------------------
LinVarChaC2AMPAS::~LinVarChaC2AMPAS() {
  mpas_linvarcha_c2a_delete_f90(keyLinVarChaC2A_);
  oops::Log::trace() << "LinVarChaC2AMPAS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2AMPAS::multiply(const IncrementMPAS & dxc,
                                IncrementMPAS & dxa) const {
  oops::Log::trace() << "LinVarChaC2AMPAS::multiply start" << std::endl;
  mpas_linvarcha_c2a_multiply_f90(keyLinVarChaC2A_, dxc.toFortran(),
                              dxa.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AMPAS::multiplyInverse(const IncrementMPAS & dxa,
                                       IncrementMPAS & dxc) const {
  oops::Log::trace() << "LinVarChaC2AMPAS::multiplyInverse start" << std::endl;
  mpas_linvarcha_c2a_multiplyinverse_f90(keyLinVarChaC2A_, dxa.toFortran(),
                              dxc.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AMPAS::multiplyAD(const IncrementMPAS & dxa,
                                  IncrementMPAS & dxc) const {
  oops::Log::trace() << "LinVarChaC2AMPAS::multiplyAD start" << std::endl;
  mpas_linvarcha_c2a_multiplyadjoint_f90(keyLinVarChaC2A_,
                                         dxa.toFortran(),
                                         dxc.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AMPAS::multiplyInverseAD(const IncrementMPAS & dxc,
                                         IncrementMPAS & dxa) const {
  oops::Log::trace() << "LinVarChaC2AMPAS::multiplyInverseAD start"
                     << std::endl;
  mpas_linvarcha_c2a_multiplyinverseadjoint_f90(keyLinVarChaC2A_,
                                                dxc.toFortran(),
                                                dxa.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AMPAS::print(std::ostream & os) const {
  os << "MPAS change variable";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
