/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_INSTANTIATEMPASVARCHANGEFACTORY_H_
#define MODEL_INSTANTIATEMPASVARCHANGEFACTORY_H_

#include "model/VarChangeMPAS.h"
#include "model/MPASTraits.h"
#include "oops/interface/LinearVariableChange.h"

namespace mpas {

void instantiateMPASVarChangeFactory() {
  static oops::LinearVariableChangeMaker<mpas::MPASTraits,
         oops::LinearVariableChange<mpas::MPASTraits, mpas::VarChangeMPAS> >
         makerVarChangeMPAS_("Control2Model");
}

}  // namespace mpas

#endif  // MODEL_INSTANTIATEMPASVARCHANGEFACTORY_H_
