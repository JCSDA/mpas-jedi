/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef MODEL_INSTANTIATELOCALIZATIONFACTORY_H_
#define MODEL_INSTANTIATELOCALIZATIONFACTORY_H_

#include "LocalizationMatrixMPAS.h"
#include "MPASTraits.h"

namespace mpas {

void instantiateLocalizationFactory() {
  // static oops::LocalizationMaker<mpas::MPASTraits, LocalizationMatrix>
  //               makerWSpeed_("MPAS");
}

}  // namespace mpas

#endif  // MODEL_INSTANTIATELOCALIZATIONFACTORY_H_

