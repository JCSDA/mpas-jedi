/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "LfricTraits.h"
#include "oops/runs/Variational.h"
#include "RunLfric.h"

int main(int argc,  char ** argv) {
  lfirc::RunLfric run(argc, argv);
  oops::Variational<lfric::LfricTraits> var;
  run.execute(var);
  return 0;
};

