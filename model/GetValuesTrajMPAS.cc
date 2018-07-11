/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>
#include "model/GetValuesTrajMPAS.h"
#include "oops/util/Logger.h"
#include "Fortran.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
GetValuesTrajMPAS::GetValuesTrajMPAS() {
  oops::Log::trace() << "GetValuesTrajMPAS constructor starting"
                     << std::endl;
  mpas_getvaltraj_setup_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTrajMPAS constructor done"
                     << keyGetValuesTraj_ << std::endl;
}
// -----------------------------------------------------------------------------
GetValuesTrajMPAS::~GetValuesTrajMPAS() {
  oops::Log::trace() << "GetValuesTrajMPAS destructor starting"
                     << std::endl;
  mpas_getvaltraj_delete_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTrajMPAS destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
