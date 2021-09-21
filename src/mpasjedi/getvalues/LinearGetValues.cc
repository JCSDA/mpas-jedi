/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/abor1_cpp.h"

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/getvalues/LinearGetValues.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/StateMPAS.h"
#include "mpasjedi/VariableChanges/Model2GeoVars/LinVarChaModel2GeoVars.h"
#include "mpasjedi/VariableChanges/Model2GeoVars/VarChaModel2GeoVars.h"

namespace mpas {

// -----------------------------------------------------------------------------

LinearGetValues::LinearGetValues(const GeometryMPAS & geom,
  const ufo::Locations & locs, const eckit::Configuration & config)
  : locs_(locs), geom_(new GeometryMPAS(geom)), linearmodel2geovars_(),
  model2geovars_() {
  oops::Log::trace() << "LinearGetValues::LinearGetValues starting"
    << std::endl;

  // Create the variable change object
  {
  util::Timer timervc(classname(), "VarChaModel2GeoVars");
  model2geovars_.reset(new VarChaModel2GeoVars(*geom_, config));
  }

  // Call GetValues consructor
  {
  util::Timer timergv(classname(), "LinearGetValues");

  mpas_lineargetvalues_create_f90(keyLinearGetValues_, geom_->toFortran(), locs,
    config);
  }

  oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

LinearGetValues::~LinearGetValues() {
  oops::Log::trace() << "LinearGetValues::~LinearGetValues starting"
    << std::endl;

  {
  util::Timer timergv(classname(), "~LinearGetValues");

  mpas_lineargetvalues_delete_f90(keyLinearGetValues_);
  }

  {
  util::Timer timervc(classname(), "~LinVarChaModel2GeoVars");

  for (lvcIter jlvc = linearmodel2geovars_.begin(); jlvc != linearmodel2geovars_.end(); ++jlvc) {
    delete jlvc->second;
  }
  }

  oops::Log::trace() << "LinearGetValues::~LinearGetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

const LinVarChaModel2GeoVars * LinearGetValues::getLinVarCha(const util::DateTime & t1) const {
  lvcIterCnst jlvc = linearmodel2geovars_.find(t1);
  if (jlvc == linearmodel2geovars_.end()) {
    oops::Log::error() << "LinearGetValues::getLinVarCha: linear variable change not available " <<
                          "at time " << t1 << std::endl;
    ABORT("LinearGetValues::getLinVarCha: linear variable change not available");
  }
  return jlvc->second;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::setTrajectory(const StateMPAS & state,
  const util::DateTime & t1, const util::DateTime & t2,
  ufo::GeoVaLs & geovals) {
  oops::Log::trace() << "LinearGetValues::setTrajectory starting" << std::endl;

  // Create state with geovals variables
  StateMPAS geovars(*geom_, geovals.getVars(), state.validTime());
  model2geovars_->changeVar(state, geovars);

  // Create the linear variable change object
  {
  util::Timer timervc(classname(), "LinearVarChaModel2GeoVars");
  ASSERT(linearmodel2geovars_.find(t1) == linearmodel2geovars_.end());
  char sep = '.';
  eckit::LocalConfiguration dummyconfig(sep);
  LinVarChaModel2GeoVars * linearmodel2geovars = new LinVarChaModel2GeoVars(state, state, *geom_,
                                                                            dummyconfig);
  linearmodel2geovars_[t1] = linearmodel2geovars;
  }

  {
  util::Timer timergv(classname(), "SetTrajectory");

  mpas_lineargetvalues_set_trajectory_f90(
    keyLinearGetValues_, geom_->toFortran(), geovars.toFortran(),
    t1, t2, locs_, geovals.toFortran());
  }

  oops::Log::trace() << "LinearGetValues::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsTL(const IncrementMPAS & inc,
  const util::DateTime & t1, const util::DateTime & t2,
  ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeoVaLsTL starting" << std::endl;

  // Create increment with geovals variables
  IncrementMPAS incgeovars(*geom_, geovals.getVars(), inc.validTime());

  {
  util::Timer timervc(classname(), "multiply");
  const LinVarChaModel2GeoVars * linearmodel2geovars = this->getLinVarCha(t1);
  linearmodel2geovars->multiply(inc, incgeovars);
  }

  {
  util::Timer timergv(classname(), "fillGeoVaLsTL");
  mpas_lineargetvalues_fill_geovals_tl_f90(
    keyLinearGetValues_, geom_->toFortran(), incgeovars.toFortran(),
    t1, t2, locs_, geovals.toFortran());
  }

  oops::Log::trace() << "LinearGetValues::fillGeoVaLsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsAD(IncrementMPAS & inc,
  const util::DateTime & t1, const util::DateTime & t2,
  const ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeoVaLsAD starting" << std::endl;

  // Create increment with geovals variables
  IncrementMPAS incgeovars(*geom_, geovals.getVars(), inc.validTime());

  {
  util::Timer timergv(classname(), "fillGeoVaLsAD");

  mpas_lineargetvalues_fill_geovals_ad_f90(
    keyLinearGetValues_, geom_->toFortran(), incgeovars.toFortran(),
    t1, t2, locs_, geovals.toFortran());
  }

  // Change variables
  {
  util::Timer timervc(classname(), "multiplyAD");
  const LinVarChaModel2GeoVars * linearmodel2geovars = this->getLinVarCha(t1);
  linearmodel2geovars->multiplyAD(incgeovars, inc);
  }

  oops::Log::trace() << "LinearGetValues::fillGeoVaLsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearGetValues::print(std::ostream & os) const {
  os << " LinearGetValues for mpas-jedi" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace mpas
