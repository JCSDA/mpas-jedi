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
#include "GetValuesTrajMPAS.h"
#include "IncrementMPAS.h"
#include "LocalizationMatrixMPAS.h"
#include "ModelBiasMPAS.h"
#include "ModelBiasIncrementMPAS.h"
#include "ModelBiasCovarianceMPAS.h"
#include "StateMPAS.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsOperator.h"
#include "ufo/LinearObsOperator.h"

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

  typedef mpas::LocalizationMatrixMPAS   LocalizationMatrix;

  typedef mpas::GetValuesTrajMPAS        InterpolatorTraj;

  typedef ufo::ObsOperator               ObsOperator;
  typedef ufo::LinearObsOperator         LinearObsOperator;
  typedef ioda::ObsSpace                 ObsSpace;
  typedef ioda::ObsVector                ObsVector;

  typedef ufo::ObsBias                   ObsAuxControl;
  typedef ufo::ObsBiasIncrement          ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance         ObsAuxCovariance;

  typedef ufo::GeoVaLs                   GeoVaLs;
  typedef ufo::Locations                 Locations;
};

}  // namespace mpas

#endif  // MODEL_MPASTRAITS_H_
