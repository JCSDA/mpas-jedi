/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/runs/Run.h>
#include <oops/runs/Variational.h>

#include <saber/oops/instantiateCovarFactory.h>
#include <saber/oops/instantiateLocalizationFactory.h>
#include <saber/oops/instantiateVariableChangeFactory.h>

#include <ufo/instantiateObsFilterFactory.h>
#include <ufo/ObsTraits.h>

#include "mpasjedi/instantiateLocalizationFactory.h"
#include "mpasjedi/MPASTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  mpas::instantiateLocalizationFactory();
  saber::instantiateCovarFactory<mpas::MPASTraits>();
  saber::instantiateLocalizationFactory<mpas::MPASTraits>();
  saber::instantiateVariableChangeFactory<mpas::MPASTraits>();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::Variational<mpas::MPASTraits, ufo::ObsTraits> var;
  return run.execute(var);
}

