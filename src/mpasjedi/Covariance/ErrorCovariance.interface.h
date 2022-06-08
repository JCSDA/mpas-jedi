/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include "mpasjedi/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mpas {

extern "C" {
  void mpas_b_setup_f90(F90bmat &, const eckit::Configuration &,
                        const F90geom &);
  void mpas_b_delete_f90(F90bmat &);

  void mpas_b_linearize_f90(const F90bmat &,
                            const eckit::Configuration &);

  void mpas_b_mult_f90(const F90bmat &, const F90inc &, const F90inc &);
  void mpas_b_invmult_f90(const F90bmat &, const F90inc &, const F90inc &);

  void mpas_b_randomize_f90(const F90bmat &, const F90inc &);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
