/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/getvalues/GetValues.h"
#include "mpasjedi/StateMPAS.h"
#include "mpasjedi/VariableChanges/Model2GeoVars/VarChaModel2GeoVars.h"

namespace mpas {

// -----------------------------------------------------------------------------

GetValues::GetValues(const GeometryMPAS & geom, const ufo::Locations & locs,
  const eckit::Configuration & config) : locs_(locs),
  geom_(new GeometryMPAS(geom)) {
  oops::Log::trace() << "GetValues::GetValues starting" << std::endl;

  // Create the variable change object
  {
  util::Timer timervc(classname(), "VarChaModel2GeoVars");
  model2geovars_.reset(new VarChaModel2GeoVars(*geom_, config));
  }

  // Call GetValues consructor
  {
  util::Timer timergv(classname(), "GetValues");

  mpas_getvalues_create_f90(keyGetValues_, geom_->toFortran(), locs_, config);
  }
  oops::Log::trace() << "GetValues::GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

GetValues::~GetValues() {
  oops::Log::trace() << "GetValues::~GetValues starting" << std::endl;
  mpas_getvalues_delete_f90(keyGetValues_);
  oops::Log::trace() << "GetValues::~GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
/*! \brief fill a GeoVaLs object using a StateMPAS object
*
* \details **fillGeoVaLs()** Calls the fortran subroutine that will convert then
* interpolate the state variable fields to the requested GeoVaLs variables then
* locations, respectively. 
*
* \param[in] state reference to the input StateMPAS object
* \param[in] t1 DateTime that is the beginning of the requested time window
* \param[in] t2 DateTime that is the end of the requested time window
* \param[out] geovals reference to the GeoVaLs object that will be populated
*
*/
void GetValues::fillGeoVaLs(const StateMPAS & state, const util::DateTime & t1,
                            const util::DateTime & t2,
                            ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "GetValues::fillGeoVaLs starting" << std::endl;

  if (geovals.getVars() <= state.variables()) {
    util::Timer timergv(classname(), "fillGeoVaLs");
    mpas_getvalues_fill_geovals_f90(keyGetValues_, geom_->toFortran(),
                                    state.toFortran(), t1, t2, locs_,
                                    geovals.toFortran());
  } else {
    // Create state with geovals variables
    StateMPAS geovars(*geom_, geovals.getVars(), state.validTime());

    {
    util::Timer timervc(classname(), "changeVar");
    model2geovars_->changeVar(state, geovars);
    }

    // Fill GeoVaLs
    mpas_getvalues_fill_geovals_f90(keyGetValues_, geom_->toFortran(),
                                    geovars.toFortran(), t1, t2, locs_,
                                    geovals.toFortran());
  }
  oops::Log::trace() << "GetValues::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

void GetValues::print(std::ostream & os) const {
  os << " GetValues for mpas-jedi" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace mpas

