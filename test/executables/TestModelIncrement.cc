/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/test/base/ModelIncrement.h>
#include "oops/runs/Run.h"

#include "mpasjedi/Traits.h"

int main(const int argc, const char ** argv) {
  oops::Run run(argc, argv);
  test::ModelIncrement<mpas::Traits> tests;
  return run.execute(tests);
}

