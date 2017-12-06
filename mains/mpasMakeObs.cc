/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/LfricTraits.h"
#include "oops/runs/MakeObs.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::MakeObs<lfric::LfricTraits> mkobs;
  run.execute(mkobs);
  return 0;
};
