/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_MPASTRAITS_H_
#define MPASJEDI_MPASTRAITS_H_

#include <string>

#include "mpasjedi/ErrorCovarianceMPAS.h"
#include "mpasjedi/GeometryMPAS.h"
#include "mpasjedi/IncrementMPAS.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/ModelBiasCovarianceMPAS.h"
#include "mpasjedi/ModelBiasIncrementMPAS.h"
#include "mpasjedi/ModelBiasMPAS.h"
#include "mpasjedi/StateMPAS.h"
#include "mpasjedi/VariableChange/VariableChange.h"

namespace mpas {

struct MPASTraits {
  static std::string name() {return "MPAS";}
  static std::string nameCovar() {return "MPASstatic";}
  static std::string nameCovar4D() {return "MPASstatic";}

  typedef mpas::ErrorCovarianceMPAS      Covariance;
  typedef mpas::GeometryMPAS             Geometry;
  typedef mpas::IncrementMPAS            Increment;
  typedef mpas::StateMPAS                State;
  typedef mpas::ModelBiasMPAS            ModelAuxControl;
  typedef mpas::ModelBiasIncrementMPAS   ModelAuxIncrement;
  typedef mpas::ModelBiasCovarianceMPAS  ModelAuxCovariance;
  typedef mpas::LinearVariableChange     LinearVariableChange;
  typedef mpas::VariableChange           VariableChange;
};

}  // namespace mpas

#endif  // MPASJEDI_MPASTRAITS_H_
