/*
 * (C) Copyright 2017-2023 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <string>

#include "mpasjedi/Covariance/ErrorCovariance.h"
#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/GeometryIterator/GeometryIterator.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/ModelBias/ModelBias.h"
#include "mpasjedi/ModelBias/ModelBiasCovariance.h"
#include "mpasjedi/ModelBias/ModelBiasIncrement.h"
#include "mpasjedi/ModelData/ModelData.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/VariableChange/VariableChange.h"

namespace mpas {

struct Traits {
  static std::string name() {return "MPAS";}
  static std::string nameCovar() {return "MPASstatic";}
  static std::string nameCovar4D() {return "MPASstatic";}

  typedef mpas::ErrorCovariance      Covariance;
  typedef mpas::Geometry             Geometry;
  typedef mpas::GeometryIterator     GeometryIterator;
  typedef mpas::Increment            Increment;
  typedef mpas::State                State;

  typedef mpas::ModelBias            ModelAuxControl;
  typedef mpas::ModelBiasIncrement   ModelAuxIncrement;
  typedef mpas::ModelBiasCovariance  ModelAuxCovariance;

  typedef mpas::LinearVariableChange LinearVariableChange;
  typedef mpas::VariableChange       VariableChange;
  typedef mpas::ModelData            ModelData;
};

}  // namespace mpas
