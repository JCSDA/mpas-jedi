/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <oops/runs/Run.h>

#include <saber/oops/instantiateCovarFactory.h>
#include <saber/oops/ProcessPerts.h>

#include "mpasjedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<mpas::Traits>();
  saber::ProcessPerts<mpas::Traits> pp;
  return run.execute(pp);
}
