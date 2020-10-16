/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/runs/HofX.h>

#include <ufo/instantiateObsFilterFactory.h>
#include <ufo/ObsTraits.h>

#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/RunMPAS.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::HofX<mpas::MPASTraits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}

