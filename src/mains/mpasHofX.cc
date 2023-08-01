/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/HofX4D.h"
#include "oops/runs/Run.h"

#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

#include "mpasjedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<mpas::Traits>();
  ufo::instantiateObsFilterFactory();
  oops::HofX4D<mpas::Traits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}

