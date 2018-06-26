/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "MPASTraits.h"
#include "instantiateLocalizationFactory.h"
#include "oops/runs/Variational.h"
#include "RunMPAS.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  mpas::instantiateLocalizationFactory();
  oops::Variational<mpas::MPASTraits> var;
  run.execute(var);
  return 0;
}

