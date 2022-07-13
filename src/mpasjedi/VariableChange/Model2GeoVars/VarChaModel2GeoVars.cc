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
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"
#include "mpasjedi/VariableChange/Model2GeoVars/VarChaModel2GeoVars.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaModel2GeoVars>
  makerVarChaModel2GeoVars_("Model2GeoVars");
static VariableChangeMaker<VarChaModel2GeoVars>
  makerVarChaDefault_("default");
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVars::VarChaModel2GeoVars(const Geometry & geom,
 const eckit::LocalConfiguration & config) : VariableChangeBase(), geom_(geom) {
  util::Timer timer(classname(), "VarChaModel2GeoVars");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_vc_model2geovars_create_f90(keyFtnConfig_, geom_.toFortran(), config);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVars::~VarChaModel2GeoVars() {
  util::Timer timer(classname(), "~VarChaModel2GeoVars");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  mpasjedi_vc_model2geovars_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVars::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar starting" << std::endl;
  mpasjedi_vc_model2geovars_changevar_f90(keyFtnConfig_, geom_.toFortran(), xin.toFortran(),
                                         xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVars::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  xout = xin;
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVars::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace mpas
