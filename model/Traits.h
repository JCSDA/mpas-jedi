/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_MPASTRAITS_H_
#define MPAS_MODEL_MPASTRAITS_H_

#include <string>

#include "model/ErrorCovariance.h"
#include "model/Geometry.h"
#include "model/Increment.h"
#include "model/State.h"
#include "model/Variables.h"

namespace mpas {

struct Traits {
  static std::string name() {return "MPAS";}
  static std::string nameCovar() {return "mpasError";}

  typedef mpas::ErrorCovariance     Covariance;
  typedef mpas::Geometry            Geometry;
  typedef mpas::Variables           Variables;
  typedef mpas::State               State;
  typedef mpas::Increment           Increment;
};

}  // namespace mpas

#endif  // MPAS_MODEL_MPASTRAITS_H_
