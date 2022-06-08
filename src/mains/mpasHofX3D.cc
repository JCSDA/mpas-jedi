/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "mpasjedi/Traits.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsErrorFactory();
  ufo::instantiateObsFilterFactory();
  oops::HofX3D<mpas::Traits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}
