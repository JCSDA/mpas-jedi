/*
 * (C) Copyright 2017-2020  UCAR.
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
#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "mpasjedi/LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaModel2GeoVars>
  makerLinVarChaModel2GeoVars_("Model2GeoVars");
static LinearVariableChangeMaker<LinVarChaModel2GeoVars>
  makerLinVarChaModel2GeoDef_("default");
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVars::LinVarChaModel2GeoVars(
  const State & bg, const State & fg, const Geometry & resol,
  const eckit::LocalConfiguration & config) : LinearVariableChangeBase(),
  geom_(resol) {
  util::Timer timer(classname(), "LinVarChaModel2GeoVars");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_lvc_model2geovars_create_f90(keyFtnConfig_, geom_.toFortran(), bg.toFortran(),
                                       fg.toFortran(), config);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVars::~LinVarChaModel2GeoVars() {
  util::Timer timer(classname(), "~LinVarChaModel2GeoVars");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  mpasjedi_lvc_model2geovars_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::changeVarTL(const Increment & dxin, Increment & dxout)
  const {
  util::Timer timer(classname(), "changeVarTL");
  oops::Log::trace() << classname() << " changeVarTL starting" << std::endl;
  mpasjedi_lvc_model2geovars_multiply_f90(keyFtnConfig_, geom_.toFortran(),
                                         dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " changeVarTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::changeVarInverseTL(const Increment & dxin, Increment & dxout)
  const {
  util::Timer timer(classname(), "changeVarInverseTL");
  oops::Log::trace() << classname() << " changeVarInverseTL starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " changeVarInverseTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::changeVarAD(const Increment & dxin, Increment & dxout)
  const {
  util::Timer timer(classname(), "changeVarAD");
  oops::Log::trace() << classname() << " changeVarAD starting" << dxin << std::endl;
  mpasjedi_lvc_model2geovars_multiplyadjoint_f90(keyFtnConfig_, geom_.toFortran(),
                                                dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " changeVarAD done" << dxout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::changeVarInverseAD(const Increment & dxin, Increment & dxout)
  const {
  util::Timer timer(classname(), "changeVarInverseAD");
  oops::Log::trace() << classname() << " changeVarInverseAD starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " changeVarInverseAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::print(std::ostream & os)
  const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace mpas
