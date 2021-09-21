/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/test/interface/LinearGetValues.h>
#include <oops/runs/Run.h>

#include <ufo/ObsTraits.h>

#include "mpasjedi/MPASTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::LinearGetValues<mpas::MPASTraits, ufo::ObsTraits> tests;
  return run.execute(tests);
}

