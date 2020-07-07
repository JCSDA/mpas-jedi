/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_MPASTRAITS_H_
#define MODEL_MPASTRAITS_H_

#include <string>

#include "ErrorCovarianceMPAS.h"
#include "GeometryMPAS.h"
#include "IncrementMPAS.h"
#include "LocalizationMatrixMPAS.h"
#include "ModelBiasMPAS.h"
#include "ModelBiasIncrementMPAS.h"
#include "ModelBiasCovarianceMPAS.h"
#include "StateMPAS.h"
#include "getvalues/GetValues.h"
#include "getvalues/LinearGetValues.h"

namespace mpas {

struct MPASTraits {
  static std::string name() {return "MPAS";}
  static std::string nameCovar() {return "MPASstatic";}
  static std::string nameCovar4D() {return "MPASstatic";}

  typedef mpas::GeometryMPAS             Geometry;

  typedef mpas::StateMPAS                State;
  typedef mpas::IncrementMPAS            Increment;
  typedef mpas::ErrorCovarianceMPAS      Covariance;

  typedef mpas::ModelBiasMPAS            ModelAuxControl;
  typedef mpas::ModelBiasIncrementMPAS   ModelAuxIncrement;
  typedef mpas::ModelBiasCovarianceMPAS  ModelAuxCovariance;

  typedef mpas::GetValues                GetValues;
  typedef mpas::LinearGetValues          LinearGetValues;

  typedef mpas::LocalizationMatrixMPAS   LocalizationMatrix;
};

}  // namespace mpas

#endif  // MODEL_MPASTRAITS_H_
