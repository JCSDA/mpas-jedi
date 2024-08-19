/*
 * (C) Copyright 2017-2023 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include "mpasjedi/Fortran.h"

// Forward declarations
namespace atlas {
  namespace field {
    class FieldSetImpl;
  }
}

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
}

namespace mpas {

extern "C" {

// -----------------------------------------------------------------------------
//  Increment
// -----------------------------------------------------------------------------
  void mpas_increment_create_f90(F90inc &, const F90geom &,
                             const oops::Variables &);
  void mpas_increment_delete_f90(F90inc &);
  void mpas_increment_copy_f90(const F90inc &, const F90inc &);
  void mpas_increment_zero_f90(const F90inc &);
  void mpas_increment_ones_f90(const F90inc &);
  void mpas_increment_axpy_inc_f90(const F90inc &, const real_type &,
                                   const F90inc &);
  void mpas_increment_axpy_state_f90(const F90inc &, const real_type &,
                                     const F90state &);
  void mpas_increment_read_file_f90(const F90inc &,
                                const eckit::Configuration &,
                                util::DateTime &);
  void mpas_increment_write_file_f90(const F90inc &,
                                 const eckit::Configuration &,
                                 const util::DateTime &);
  void mpas_increment_gpnorm_f90(const F90inc &, const int &, real_type &);
  void mpas_increment_rms_f90(const F90inc &, real_type &);
  void mpas_increment_diff_incr_f90(const F90inc &, const F90state &,
                                    const F90state &);
  void mpas_increment_self_add_f90(const F90inc &, const F90inc &);
  void mpas_increment_self_sub_f90(const F90inc &, const F90inc &);
  void mpas_increment_self_mul_f90(const F90inc &, const real_type &);
  void mpas_increment_dot_prod_f90(const F90inc &, const F90inc &, real_type &);
  void mpas_increment_self_schur_f90(const F90inc &, const F90inc &);
  void mpas_increment_random_f90(const F90inc &);
  void mpas_increment_to_fieldset_f90(const F90inc &,
                                   const F90geom &,
                                   const oops::Variables &,
                                   atlas::field::FieldSetImpl *,
                                   const bool &,
                                   const bool &);
  void mpas_increment_from_fieldset_f90(const F90inc &,
                                     const F90geom &,
                                     const oops::Variables &,
                                     const atlas::field::FieldSetImpl *,
                                     const bool &,
                                     const bool &);
  void mpas_increment_dirac_f90(const F90inc &,
                                const eckit::Configuration &);
  void mpas_increment_sizes_f90(const F90inc &, int &, int &);
  void mpas_increment_serial_size_f90(const F90inc &, std::size_t &);
  void mpas_increment_serialize_f90(const F90inc &, const std::size_t &,
                                    real_type[]);
  void mpas_increment_deserialize_f90(const F90inc &, const std::size_t &,
                                      const real_type[], const std::size_t &);
  void mpas_increment_getpoint_f90(const F90inc &, const F90iter &, double &, const int &);
  void mpas_increment_setpoint_f90(F90inc &, const F90iter &, const double &, const int &);
};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace mpas
