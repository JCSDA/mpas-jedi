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

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/LinearVariableChange/Control2Analysis/LinVarChaC2A.h"
#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/StateMPAS.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaC2A>
  makerLinVarChaC2A_("Control2Analysis");
// -----------------------------------------------------------------------------
LinVarChaC2A::LinVarChaC2A(const StateMPAS & bg,
                           const StateMPAS & fg,
                           const GeometryMPAS & resol,
                           const eckit::LocalConfiguration & config)
: LinearVariableChangeBase(), geom_(new GeometryMPAS(resol)), vars_(config, "input variables") {
  util::Timer timer(classname(), "LinVarChaC2A");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_->toFortran(),
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
void LinVarChaC2A::multiply(const IncrementMPAS & dxc, IncrementMPAS & dxa)
  const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                     dxc.toFortran(), dxa.toFortran());
  oops::Log::trace() << classname() << " multiply done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverse(const IncrementMPAS & dxa, IncrementMPAS & dxc)
  const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
  oops::Log::trace() << classname() << " multiplyInverse done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::multiplyAD(const IncrementMPAS & dxa, IncrementMPAS & dxc)
  const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
  oops::Log::trace() << classname() << " multiplyAD done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverseAD(const IncrementMPAS & dxc, IncrementMPAS & dxa)
  const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  mpasjedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                   dxc.toFortran(), dxa.toFortran());
  oops::Log::trace() << classname() << " multiplyInverseAD done" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2A::print(std::ostream & os)
  const {
  os << classname() << " variable change";
}
// -----------------------------------------------------------------------------
}  // namespace mpas
