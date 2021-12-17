/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mpasjedi/Fortran.h"

namespace mpas {
  typedef int F90lvc_M2G;
  extern "C" {
  void mpasjedi_lvc_model2geovars_create_f90(const F90lvc_M2G &, const F90geom &, const F90state &,
                                            const F90state &, const eckit::LocalConfiguration &);
  void mpasjedi_lvc_model2geovars_delete_f90(F90lvc_M2G &);
  void mpasjedi_lvc_model2geovars_multiply_f90(const F90lvc_M2G &, const F90geom &, const F90inc &,
                                              const F90inc &);
  void mpasjedi_lvc_model2geovars_multiplyadjoint_f90(const F90lvc_M2G &, const F90geom &,
                                                     const F90inc &, const F90inc &);
  }  // extern "C"
}  // namespace mpas
