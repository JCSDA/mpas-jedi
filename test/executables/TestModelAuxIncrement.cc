/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/test/interface/ModelAuxIncrement.h>

#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/RunMPAS.h"

int main(const int argc, const char ** argv) {
  mpas::RunMPAS run(argc, argv);
  test::ModelAuxIncrement<mpas::MPASTraits> tests;
  return run.execute(tests);
}

