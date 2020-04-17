/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "RunMPAS.h"
#include "oops/runs/EnsHofX.h"
#include "MPASTraits.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  oops::EnsHofX<mpas::MPASTraits> enshofx;
  return run.execute(enshofx);
}

