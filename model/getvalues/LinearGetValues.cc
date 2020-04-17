/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/getvalues/LinearGetValues.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

LinearGetValues::LinearGetValues(const GeometryMPAS & geom, const ufo::Locations & locs)  : locs_(locs),
     geom_(new GeometryMPAS(geom)) {
  oops::Log::trace() << "LinearGetValues::LinearGetValues starting" << std::endl;
  mpas_lineargetvalues_create_f90(keyLinearGetValues_, geom.toFortran(), locs.toFortran());
  oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

LinearGetValues::~LinearGetValues() {
  oops::Log::trace() << "LinearGetValues::~LinearGetValues starting" << std::endl;

  mpas_lineargetvalues_delete_f90(keyLinearGetValues_);

  oops::Log::trace() << "LinearGetValues::~LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::setTrajectory(const StateMPAS & state, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) {
  oops::Log::trace() << "LinearGetValues::setTrajectory starting" << std::endl;

  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;

  mpas_lineargetvalues_set_trajectory_f90(keyLinearGetValues_, geom_->toFortran(),
                                          state.toFortran(), &t1p, &t2p,
                                          locs_.toFortran(), geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsTL(const IncrementMPAS & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeoVaLsTL starting" << std::endl;

  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  mpas_lineargetvalues_fill_geovals_tl_f90(keyLinearGetValues_, geom_->toFortran(),
                                           inc.toFortran(), &t1p, &t2p,
                                           locs_.toFortran(), geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::fillGeoVaLsTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsAD(IncrementMPAS & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, const ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeoVaLsAD starting" << std::endl;

  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  mpas_lineargetvalues_fill_geovals_ad_f90(keyLinearGetValues_, geom_->toFortran(),
                                           inc.toFortran(), &t1p, &t2p,
                                           locs_.toFortran(), geovals.toFortran());

  oops::Log::trace() << "LinearGetValues::fillGeoVaLsAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::print(std::ostream & os) const {
  os << " LinearGetValues for mpas-jedi" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
