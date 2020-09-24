/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_INCREMENTMPASFORTRAN_H_
#define MODEL_INCREMENTMPASFORTRAN_H_

#include "Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace util {
  class DateTime;
  class Duration;
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
  void mpas_increment_axpy_inc_f90(const F90inc &, const double &,
                                   const F90inc &);
  void mpas_increment_axpy_state_f90(const F90inc &, const double &,
                                     const F90state &);

  void mpas_increment_change_resol_f90(const F90inc &, const F90inc &);
  void mpas_increment_read_file_f90(const F90inc &,
                                const eckit::Configuration &,
                                util::DateTime &);
  void mpas_increment_write_file_f90(const F90inc &,
                                 const eckit::Configuration &,
                                 const util::DateTime &);
  void mpas_increment_gpnorm_f90(const F90inc &, const int &, double &);
  void mpas_increment_rms_f90(const F90inc &, double &);
  void mpas_increment_diff_incr_f90(const F90inc &, const F90state &,
                                    const F90state &);
  void mpas_increment_self_add_f90(const F90inc &, const F90inc &);
  void mpas_increment_self_sub_f90(const F90inc &, const F90inc &);
  void mpas_increment_self_mul_f90(const F90inc &, const double &);
  void mpas_increment_dot_prod_f90(const F90inc &, const F90inc &, double &);
  void mpas_increment_self_schur_f90(const F90inc &, const F90inc &);
  void mpas_increment_random_f90(const F90inc &);
  void mpas_increment_set_atlas_f90(const F90inc &,
                                    const F90geom &,
                                    const oops::Variables &,
                                    const util::DateTime * const *,
                                    atlas::field::FieldSetImpl *);
  void mpas_increment_to_atlas_f90(const F90inc &,
                                   const F90geom &,
                                   const oops::Variables &,
                                   const util::DateTime * const *,
                                   atlas::field::FieldSetImpl *);
  void mpas_increment_from_atlas_f90(const F90inc &,
                                     const F90geom &,
                                     const oops::Variables &,
                                     const util::DateTime * const *,
                                     atlas::field::FieldSetImpl *);
  void mpas_increment_dirac_f90(const F90inc &,
                                const eckit::Configuration &);
  void mpas_increment_sizes_f90(const F90inc &, int &, int &);
  void mpas_increment_serial_size_f90(const F90inc &, int &);
  void mpas_increment_serialize_f90(const F90inc &, const std::size_t &, double[]);
  void mpas_increment_deserialize_f90(const F90inc &, const std::size_t &, const double[],
                                      const std::size_t &);

};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_INCREMENTMPASFORTRAN_H_
