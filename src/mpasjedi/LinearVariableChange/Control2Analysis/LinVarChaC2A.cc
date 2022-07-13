/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/LinearVariableChange/Control2Analysis/LinVarChaC2A.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaC2A>
  makerLinVarChaC2A_("Control2Analysis");
// -----------------------------------------------------------------------------
LinVarChaC2A::LinVarChaC2A(const State & bg,
                           const State & fg,
                           const Geometry & resol,
                           const eckit::LocalConfiguration & config)
: LinearVariableChangeBase(), geom_(resol), vars_(config, "input variables") {
  util::Timer timer(classname(), "LinVarChaC2A");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_.toFortran(),
                                   bg.toFortran(), fg.toFortran(),
                                   config, vars_);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -----------------------------------------------------------------------------
LinVarChaC2A::~LinVarChaC2A() {
  util::Timer timer(classname(), "~LinVarChaC2A");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  mpasjedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::changeVarTL(const Increment & dxc, Increment & dxa)
  const {
  util::Timer timer(classname(), "changeVarTL");
  oops::Log::trace() << classname() << " changeVarTL starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiply_f90(keyFtnConfig_, geom_.toFortran(),
                                     dxc.toFortran(), dxa.toFortran());
  oops::Log::trace() << classname() << " changeVarTL done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::changeVarInverseTL(const Increment & dxa, Increment & dxc)
  const {
  util::Timer timer(classname(), "changeVarInverseTL");
  oops::Log::trace() << classname() << " changeVarInverseTL starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_.toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
  oops::Log::trace() << classname() << " changeVarInverseTL done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::changeVarAD(const Increment & dxa, Increment & dxc)
  const {
  util::Timer timer(classname(), "changeVarAD");
  oops::Log::trace() << classname() << " changeVarAD starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyadjoint_f90(keyFtnConfig_, geom_.toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
  oops::Log::trace() << classname() << " changeVarAD done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::changeVarInverseAD(const Increment & dxc, Increment & dxa)
  const {
  util::Timer timer(classname(), "changeVarInverseAD");
  oops::Log::trace() << classname() << " changeVarInverseAD starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_.toFortran(),
                                                   dxc.toFortran(), dxa.toFortran());
  oops::Log::trace() << classname() << " changeVarInverseAD done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::print(std::ostream & os)
  const {
  os << classname() << " variable change";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
