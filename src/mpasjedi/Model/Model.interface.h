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
  void mpas_model_setup_f90(const eckit::Configuration &,
                            const F90geom &, F90model &);
  void mpas_model_delete_f90(F90model &);

  void mpas_model_prepare_integration_f90(const F90model &, const F90state &);
  void mpas_model_prepare_integration_tl_f90(const F90model &, const F90inc &);
  void mpas_model_prepare_integration_ad_f90(const F90model &, const F90inc &);

  void mpas_model_propagate_f90(const F90model &, const F90state &);
  void mpas_model_prop_traj_f90(const F90model &, const F90state &, F90traj &);
  void mpas_model_propagate_tl_f90(const F90model &, const F90inc &,
                                   const F90traj &);
  void mpas_model_propagate_ad_f90(const F90model &, const F90inc &,
                                   const F90traj &);

  void mpas_model_wipe_traj_f90(F90traj &);
  void mpas_traj_minmaxrms_f90(const F90traj &, double &);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
