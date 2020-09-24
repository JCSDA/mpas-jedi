/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_STATEMPASFORTRAN_H_
#define MODEL_STATEMPASFORTRAN_H_

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
//  State
// -----------------------------------------------------------------------------
  void mpas_state_create_f90(F90state &, const F90geom &,
                             const oops::Variables &,
                             const oops::Variables &);
  void mpas_state_delete_f90(F90state &);
  void mpas_state_copy_f90(const F90state &, const F90state &);
  void mpas_state_zero_f90(const F90state &);
  void mpas_state_axpy_f90(const F90state &, const double &, const F90state &);
  void mpas_state_add_incr_f90(const F90state &, const F90inc &);
  void mpas_state_change_resol_f90(const F90state &, const F90state &);
  void mpas_state_serial_size_f90(const F90state &, int &);
  void mpas_state_serialize_f90(const F90state &, const std::size_t &, double[]);
  void mpas_state_deserialize_f90(const F90state &, const std::size_t &, const double[],
                                      const std::size_t &);
  void mpas_state_read_file_f90(const F90state &,
                                const eckit::Configuration &,
                                util::DateTime &);
  void mpas_state_write_file_f90(const F90state &,
                                 const eckit::Configuration &,
                                 const util::DateTime &);
  void mpas_state_gpnorm_f90(const F90state &, const int &, double &);
  void mpas_state_rms_f90(const F90state &, double &);
  void mpas_state_analytic_init_f90(const F90state &, const F90geom &,
                                        const eckit::Configuration &,
                                        util::DateTime &);
  void mpas_state_sizes_f90(const F90state &, int &, int &);


};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_STATEMPASFORTRAN_H_
