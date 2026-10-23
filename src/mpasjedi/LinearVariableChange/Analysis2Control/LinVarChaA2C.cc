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
#include "mpasjedi/LinearVariableChange/Analysis2Control/LinVarChaA2C.h"
#include "mpasjedi/State/State.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaA2C>
  makerLinVarChaA2C_("Analysis2Control");
// -----------------------------------------------------------------------------
LinVarChaA2C::LinVarChaA2C(const State & bg,
                           const State & fg,
                           const Geometry & resol,
                           const eckit::LocalConfiguration & config)
: LinearVariableChangeBase(), geom_(resol), vars_(config, "input variables") {
  util::Timer timer(classname(), "LinVarChaA2C");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_.toFortran(),
                                   bg.toFortran(), fg.toFortran(),
                                   config, vars_);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -----------------------------------------------------------------------------
LinVarChaA2C::~LinVarChaA2C() {
  util::Timer timer(classname(), "~LinVarChaA2C");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  mpasjedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaA2C::changeVarTL(const Increment & dxa, Increment & dxc)
  const {
  util::Timer timer(classname(), "changeVarTL");
  oops::Log::trace() << classname() << " changeVarTL starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_.toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
  oops::Log::trace() << classname() << " changeVarInverseTL done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaA2C::changeVarAD(const Increment & dxc, Increment & dxa)
  const {
  util::Timer timer(classname(), "changeVarAD");
  oops::Log::trace() << classname() << " changeVarAD starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_.toFortran(),
                                                   dxc.toFortran(), dxa.toFortran());
  oops::Log::trace() << classname() << " changeVarInverseAD done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaA2C::print(std::ostream & os)
  const {
  os << classname() << " variable change";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
