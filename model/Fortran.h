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

// Geometry key type
typedef int F90geom;
// Model key type
typedef int F90model;
// Variables key type
typedef int F90vars;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Fields key type
typedef int F90flds;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// Observation vector key type
typedef int F90ovec;
// Obs operator key type
typedef int F90hop;
// Observation data base type
typedef int F90odb;
// Localization matrix
typedef int F90lclz;

/// Interface to Fortran MPAS model
/*!
 * The core of the MPAS model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mpas_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
  void mpas_geo_clone_f90(const F90geom &, F90geom &);
  void mpas_geo_info_f90(const F90geom &, int &, int &, int &, int &, int &, int &, int &, int &);
  void mpas_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void mpas_model_setup_f90(const eckit::Configuration * const *, const F90geom &, F90model &);
  void mpas_model_delete_f90(F90model &);

  void mpas_model_prepare_integration_f90(const F90model &, const F90flds &);
  void mpas_model_prepare_integration_tl_f90(const F90model &, const F90flds &);
  void mpas_model_prepare_integration_ad_f90(const F90model &, const F90flds &);

  void mpas_model_propagate_f90(const F90model &, const F90flds &);
  void mpas_model_prop_traj_f90(const F90model &, const F90flds &, F90traj &);
  void mpas_model_propagate_tl_f90(const F90model &, const F90flds &, const F90traj &);
  void mpas_model_propagate_ad_f90(const F90model &, const F90flds &, const F90traj &);

  void mpas_model_wipe_traj_f90(F90traj &);
  void mpas_traj_minmaxrms_f90(const F90traj &, double &);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void mpas_field_create_f90(F90flds &, const F90geom &, const eckit::Configuration * const *);
  void mpas_field_delete_f90(F90flds &);

  void mpas_field_copy_f90(const F90flds &, const F90flds &);
  void mpas_field_zero_f90(const F90flds &);
  void mpas_field_self_add_f90(const F90flds &, const F90flds &);
  void mpas_field_self_sub_f90(const F90flds &, const F90flds &);
  void mpas_field_self_mul_f90(const F90flds &, const double &);
  void mpas_field_axpy_f90(const F90flds &, const double &, const F90flds &);
  void mpas_field_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void mpas_field_self_schur_f90(const F90flds &, const F90flds &);
  void mpas_field_random_f90(const F90flds &);

  void mpas_field_add_incr_f90(const F90flds &, const F90flds &);
  void mpas_field_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);

  void mpas_field_change_resol_f90(const F90flds &, const F90flds &);

  void mpas_field_read_file_f90(const F90flds &, const eckit::Configuration * const *,
                              util::DateTime * const *);
  void mpas_field_write_file_f90(const F90flds &, const eckit::Configuration * const *,
                               const util::DateTime * const *);

  void mpas_field_interp_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *, const F90goms &);
  void mpas_field_interp_tl_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *, const F90goms &);
  void mpas_field_interp_ad_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *, const F90goms &);
  void mpas_field_define_f90(const F90flds &, const int &);
  void mpas_field_convert_to_f90(const F90flds &, const int &);
  void mpas_field_convert_from_f90(const F90flds &, const int &);

  void mpas_field_gpnorm_f90(const F90flds &, const int &, double &);
  void mpas_field_sizes_f90(const F90flds &, int &, int &);
  void mpas_field_rms_f90(const F90flds &, double &);

  void mpas_field_dirac_f90(const F90flds &, const eckit::Configuration * const *);

// -----------------------------------------------------------------------------
//  Background error
// -----------------------------------------------------------------------------
  void mpas_b_setup_f90(F90bmat &, const eckit::Configuration * const *, const F90geom &);
  void mpas_b_delete_f90(F90bmat &);

  void mpas_b_linearize_f90(const F90bmat &, const eckit::Configuration * const *);

  void mpas_b_mult_f90(const F90bmat &, const F90flds &, const F90flds &);
  void mpas_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);

  void mpas_b_randomize_f90(const F90bmat &, const F90flds &);

// -----------------------------------------------------------------------------
//  Localization matrix
// -----------------------------------------------------------------------------
  void mpas_localization_setup_f90(F90lclz &, const eckit::Configuration * const *,
                                 const F90geom &);
  void mpas_localization_delete_f90(F90lclz &);
  void mpas_localization_mult_f90(const F90lclz &, const F90flds &);

// -----------------------------------------------------------------------------
//  Run
// -----------------------------------------------------------------------------
  void mpas_init_f90(const eckit::Configuration * const *);
  void mpi_finalize_f90();

// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void mpas_loc_delete_f90(F90locs &);
  void mpas_loc_nobs_f90(const F90locs &, int &);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void mpas_gom_create_f90(F90goms &);
  void mpas_gom_delete_f90(F90goms &);
  void mpas_gom_zero_f90(const F90goms &);
  void mpas_gom_random_f90(const F90goms &);
  void mpas_gom_dotprod_f90(const F90goms &, const F90goms &, double &);
  void mpas_gom_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void mpas_gom_read_file_f90(const F90goms &, const eckit::Configuration * const *);
  void mpas_gom_write_file_f90(const F90goms &, const eckit::Configuration * const *);

// -----------------------------------------------------------------------------
//  Streamfunction observations
// -----------------------------------------------------------------------------
  void mpas_stream_setup_f90(F90hop &, const eckit::Configuration * const *);
  void mpas_stream_delete_f90(F90hop &);

  void mpas_stream_equiv_f90(const F90goms &, const F90ovec &, const double &);
  void mpas_stream_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void mpas_stream_equiv_ad_f90(const F90goms &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind observations
// -----------------------------------------------------------------------------
  void mpas_wind_setup_f90(F90hop &, const eckit::Configuration * const *);
  void mpas_wind_delete_f90(F90hop &);

  void mpas_wind_equiv_f90(const F90goms &, F90ovec &, const double &);
  void mpas_wind_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void mpas_wind_equiv_ad_f90(const F90goms &, const F90ovec &, double &);

// -----------------------------------------------------------------------------
//  Wind speed observations
// -----------------------------------------------------------------------------
  void mpas_wspeed_setup_f90(F90hop &, const eckit::Configuration * const *);
  void mpas_wspeed_delete_f90(F90hop &);

  void mpas_wspeed_eqv_f90(const F90goms &, const F90ovec &, const double &);
  void mpas_wspeed_equiv_tl_f90(const F90goms &, const F90ovec &, const F90goms &, const double &);
  void mpas_wspeed_equiv_ad_f90(const F90goms &, const F90ovec &, const F90goms &, double &);

  void mpas_wspeed_gettraj_f90(const F90hop &, const int &, F90goms &);
  void mpas_wspeed_settraj_f90(const F90goms &, const F90goms &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void mpas_obsvec_setup_f90(F90ovec &, const int &, const int &);
  void mpas_obsvec_clone_f90(const F90ovec &, F90ovec &);
  void mpas_obsvec_delete_f90(F90ovec &);

  void mpas_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void mpas_obsvec_zero_f90(const F90ovec &);
  void mpas_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void mpas_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void mpas_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void mpas_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void mpas_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void mpas_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void mpas_obsvec_invert_f90(const F90ovec &);
  void mpas_obsvec_random_f90(const F90ovec &);
  void mpas_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void mpas_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void mpas_obsvec_nobs_f90(const F90ovec &, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void mpas_obsdb_setup_f90(F90odb &, const eckit::Configuration * const *);
  void mpas_obsdb_delete_f90(F90odb &);
  void mpas_obsdb_get_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void mpas_obsdb_put_f90(const F90odb &, const int &, const char *,
                        const int &, const char *, const F90ovec &);
  void mpas_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              F90locs &);
  void mpas_obsdb_getgom_f90(const F90odb &, const int &, const char *, const F90vars *,
                           const util::DateTime * const *, const util::DateTime * const *,
                           F90goms &);
  void mpas_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void mpas_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);
  void mpas_obsoper_inputs_f90(const F90hop &, F90vars *);
}

// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_MODEL_MPASFORTRAN_H_
