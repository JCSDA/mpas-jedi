/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "mpasjedi/Traits.h"

#include "oops/base/instantiateInflationFactory.h"
#include "oops/runs/EnsembleInflation.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::instantiateInflationFactory<mpas::Traits>();
  oops::Run run(argc, argv);
  oops::EnsembleInflation<mpas::Traits> ensinfl;
  return run.execute(ensinfl);
}
