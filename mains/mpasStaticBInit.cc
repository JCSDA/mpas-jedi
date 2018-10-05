/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "RunMPAS.h"
#include "oops/runs/StaticBInit.h"
#include "MPASTraits.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  oops::StaticBInit<mpas::MPASTraits> bmat;
  run.execute(bmat);
  return 0;
}
