/*
 * (C) Copyright 2017 UCAR
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
//  State
// -----------------------------------------------------------------------------
  void mpas_state_create_f90(F90state &, const F90geom &,
                             const oops::Variables &,
                             const oops::Variables &);
  void mpas_state_delete_f90(F90state &);
  void mpas_state_copy_f90(const F90state &, const F90state &);
  void mpas_state_zero_f90(const F90state &);
  void mpas_state_axpy_f90(const F90state &, const real_type &, const F90state &);
  void mpas_state_add_incr_f90(const F90state &, const F90inc &);
  void mpas_state_serial_size_f90(const F90state &, std::size_t &);
  void mpas_state_serialize_f90(const F90state &, const std::size_t &,
                                real_type[]);
  void mpas_state_deserialize_f90(const F90state &, const std::size_t &,
                                  const real_type[], const std::size_t &);
  void mpas_state_read_file_f90(const F90state &,
                                const eckit::Configuration &,
                                util::DateTime &);
  void mpas_state_write_file_f90(const F90state &,
                                 const eckit::Configuration &,
                                 const util::DateTime &);
  void mpas_state_gpnorm_f90(const F90state &, const int &, real_type &);
  void mpas_state_rms_f90(const F90state &, real_type &);
  void mpas_state_analytic_init_f90(const F90state &,
                                    const eckit::Configuration &,
                                    util::DateTime &);
  void mpas_state_sizes_f90(const F90state &, int &, int &);
  void mpas_state_to_fieldset_f90(const F90inc &, const F90geom &, const oops::Variables &,
                               atlas::field::FieldSetImpl *, const bool &, const bool &);
  void mpas_state_from_fieldset_f90(const F90inc &, const F90geom &, const oops::Variables &,
                               const atlas::field::FieldSetImpl *, const bool &, const bool &);

};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace mpas
