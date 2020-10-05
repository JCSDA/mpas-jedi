/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "RunMPAS.h"
#include "oops/runs/EnsembleApplication.h"
#include "oops/runs/HofX.h"
#include "MPASTraits.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  mpas::RunMPAS run(argc, argv);
  oops::EnsembleApplication<oops::HofX <mpas::MPASTraits, ufo::ObsTraits> >enshofx;
  return run.execute(enshofx);
}
