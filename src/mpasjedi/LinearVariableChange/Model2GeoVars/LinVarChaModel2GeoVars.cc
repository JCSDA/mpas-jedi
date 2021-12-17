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

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "mpasjedi/LinearVariableChange/Model2GeoVars/LinVarChaModel2GeoVars.h"
#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/StateMPAS.h"

namespace mpas {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaModel2GeoVars>
  makerLinVarChaModel2GeoVars_("Model2GeoVars");
static LinearVariableChangeMaker<LinVarChaModel2GeoVars>
  makerLinVarChaModel2GeoDef_("default");
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVars::LinVarChaModel2GeoVars(
  const StateMPAS & bg, const StateMPAS & fg, const GeometryMPAS & resol,
  const eckit::LocalConfiguration & config) : LinearVariableChangeBase(),
  geom_(new GeometryMPAS(resol)) {
  util::Timer timer(classname(), "LinVarChaModel2GeoVars");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  mpasjedi_lvc_model2geovars_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
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
void LinVarChaModel2GeoVars::multiply(const IncrementMPAS & dxin, IncrementMPAS & dxout)
  const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  mpasjedi_lvc_model2geovars_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                         dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiply done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::multiplyInverse(const IncrementMPAS & dxin, IncrementMPAS & dxout)
  const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " multiplyInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::multiplyAD(const IncrementMPAS & dxin, IncrementMPAS & dxout)
  const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << dxin << std::endl;
  mpasjedi_lvc_model2geovars_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD done" << dxout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::multiplyInverseAD(const IncrementMPAS & dxin, IncrementMPAS & dxout)
  const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " multiplyInverseAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVars::print(std::ostream & os)
  const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace mpas
