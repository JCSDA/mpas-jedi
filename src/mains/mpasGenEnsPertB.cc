/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/runs/GenEnsPertB.h>
#include <oops/runs/Run.h>
#include "saber/oops/instantiateCovarFactory.h"

#include "mpasjedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<mpas::Traits>();
  oops::GenEnsPertB<mpas::Traits> ensgen;
  return run.execute(ensgen);
}
