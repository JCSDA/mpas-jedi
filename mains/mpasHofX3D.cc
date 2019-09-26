/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "MPASTraits.h"
#include "RunMPAS.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "oops/runs/HofX3D.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  ufo::instantiateObsFilterFactory<mpas::MPASTraits>();
  oops::HofX3D<mpas::MPASTraits> hofx3d;
  run.execute(hofx3d);
  return 0;
}

