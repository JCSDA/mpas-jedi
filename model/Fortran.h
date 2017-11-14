/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_MPASFORTRAN_H_
#define MPAS_MODEL_MPASFORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace mpas {

/// Interface to Fortran MPAS model
/*!
 * The core of the MPAS model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mpas_geo_setup_f90(int & keyGeom, const eckit::Configuration * const *);
  void mpas_geo_clone_f90(const int & keyGeom, int & keyGeom_other);
  void mpas_geo_info_f90(const int & keyGeom, int &, int &, int &, int &);
  void mpas_geo_delete_f90(int & keyGeom);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void mpas_field_create_f90(int & keyFlds, const int &, const int & keyVars);
  void mpas_field_delete_f90(int & keyFlds);

  void mpas_field_copy_f90(const int & keyFlds, const int & keyFldsOther);
  void mpas_field_zero_f90(const int & keyFlds);
  void mpas_field_dirac_f90(const int & keyFlds, const eckit::Configuration * const *);
  void mpas_field_self_add_f90(const int & keyFlds, const int & keyFldsOther);
  void mpas_field_self_sub_f90(const int & keyFlds, const int & keyFldsOther);
  void mpas_field_self_mul_f90(const int & keyFlds, const double &);
  void mpas_field_axpy_f90(const int & keyFlds, const double &, const int & keyFldsOther);
  void mpas_field_dot_prod_f90(const int & keyFlds, const int & keyFldsOther, double &);
  void mpas_field_self_schur_f90(const int & keyFlds, const int & keyFldsOther);
  void mpas_field_random_f90(const int & keyFlds);

  void mpas_field_add_incr_f90(const int & keyFlds, const int & keyFldsOther);
  void mpas_field_diff_incr_f90(const int & keyFlds, const int & keyFldsOther,
                              const int & keyFldsOther2);

  void mpas_field_change_resol_f90(const int & keyFlds, const int & keyFldsOther);

  void mpas_field_read_file_f90(const int & keyFlds, const eckit::Configuration * const *,
                              util::DateTime * const *);
  void mpas_field_write_file_f90(const int & keyFlds, const eckit::Configuration * const *,
                               const util::DateTime * const *);

  void mpas_field_convert_to_f90(const int &, const int &);
  void mpas_field_convert_from_f90(const int &, const int &);

  void mpas_field_gpnorm_f90(const int & keyFlds, const int &, double &);
  void mpas_field_sizes_f90(const int & keyFlds, int &, int &, int &);
  void mpas_field_rms_f90(const int & keyFlds, double &);

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void mpas_var_create_f90(int & keyVars, const eckit::Configuration * const *);
  void mpas_var_clone_f90(const int & keyVars, int & keyVars_other);
  void mpas_var_info_f90(const int & keyVars, int &, int &);
  void mpas_var_delete_f90(int & keyVars);
}

// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_MPASFORTRAN_H_
