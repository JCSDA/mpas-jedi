/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/runs/HofX3D.h>
#include <oops/runs/Run.h>

#include <ufo/instantiateObsFilterFactory.h>
#include <ufo/ObsTraits.h>

#include "mpasjedi/MPASTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::HofX3D<mpas::MPASTraits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}

