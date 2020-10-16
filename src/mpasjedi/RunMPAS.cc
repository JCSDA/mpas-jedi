/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "mpasjedi/RunMPAS.h"

#include <fstream>

#include <eckit/config/Configuration.h>

#include <oops/util/Logger.h>

#include "mpasjedi/Fortran.h"

namespace mpas {

// -----------------------------------------------------------------------------

RunMPAS::RunMPAS(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating RunMPAS" << std::endl;
  mpas_run_init_f90(config());
  oops::Log::trace() << "RunMPAS created" << std::endl;
}

// -----------------------------------------------------------------------------

RunMPAS::~RunMPAS() {
  oops::Log::trace() << "Destructing RunMPAS" << std::endl;
  mpas_run_final_f90();
  oops::Log::trace() << "RunMPAS: MPI finalized" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace mpas
