/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mpasjedi/Fortran.h"

namespace mpas {

  extern "C" {
    void mpas_geom_iter_setup_f90(F90iter &, const F90geom &,
                                  const int &, const int &);
    void mpas_geom_iter_clone_f90(F90iter &, const F90iter &);
    void mpas_geom_iter_delete_f90(F90iter &);
    void mpas_geom_iter_equals_f90(const F90iter &, const F90iter&, bool &);
    void mpas_geom_iter_current_f90(const F90iter &,
                                    double &, double &, double &);
    void mpas_geom_iter_next_f90(const F90iter &);
    void mpas_geom_iter_dimension_f90(const F90iter &, int &);
  }
}  // namespace mpas
