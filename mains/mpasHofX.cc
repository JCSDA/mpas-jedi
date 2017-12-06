/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "RunLfric.h"
#include "oops/runs/HofX.h"
#include "LfricTraits.h"

int main(int argc,  char ** argv) {
  lfric::RunLfric run(argc, argv);
  oops::HofX<lfric::LfricTraits> hofx;
  run.execute(hofx);
  return 0;
};

