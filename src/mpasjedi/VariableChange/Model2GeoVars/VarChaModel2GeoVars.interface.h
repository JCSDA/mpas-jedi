/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mpasjedi/Fortran.h"

namespace mpas {
  typedef int F90vc_M2G;
  extern "C" {
  void mpasjedi_vc_model2geovars_create_f90(const F90vc_M2G &, const F90geom &,
                                            const eckit::LocalConfiguration &);
  void mpasjedi_vc_model2geovars_delete_f90(F90vc_M2G &);
  void mpasjedi_vc_model2geovars_changevar_f90(const F90vc_M2G &, const F90geom &, const F90state &,
                                              const F90state &);
  }  // extern "C"
}  // namespace mpas
