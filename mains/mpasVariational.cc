/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/runs/Variational.h"
#include "ufo/instantiateObsFilterFactory.h"

#include "MPASTraits.h"
#include "instantiateLocalizationFactory.h"
#include "instantiateMPASVarChangeFactory.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"
#include "saber/oops/instantiateVariableChangeFactory.h"
#include "RunMPAS.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  mpas::instantiateLocalizationFactory();
  mpas::instantiateMPASVarChangeFactory();
  saber::instantiateCovarFactory<mpas::MPASTraits>();
  saber::instantiateLocalizationFactory<mpas::MPASTraits>();
  saber::instantiateVariableChangeFactory<mpas::MPASTraits>();
  ufo::instantiateObsFilterFactory<mpas::MPASTraits>();
  oops::Variational<mpas::MPASTraits> var;
  run.execute(var);
  return 0;
}

