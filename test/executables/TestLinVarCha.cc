/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <oops/test/interface/LinearVariableChange.h>

#include <saber/oops/instantiateVariableChangeFactory.h>

#include "mpasjedi/instantiateMPASVarChangeFactory.h"
#include "mpasjedi/MPASTraits.h"
#include "mpasjedi/RunMPAS.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  mpas::instantiateMPASVarChangeFactory();
  saber::instantiateVariableChangeFactory<mpas::MPASTraits>();
  test::LinearVariableChange<mpas::MPASTraits> tests;
  return run.execute(tests);
}

