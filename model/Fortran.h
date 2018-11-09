/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MODEL_FORTRAN_H_
#define MODEL_FORTRAN_H_

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
// ObOp trajectory
typedef int F90ootrj;
// VarChange key
typedef int F90vc;

/// Interface to Fortran MPAS model
/*!
 * The core of the MPAS model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Run
// -----------------------------------------------------------------------------
  void mpas_run_init_f90(const eckit::Configuration * const *);
  void mpas_run_final_f90();

// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mpas_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
  void mpas_geo_clone_f90(const F90geom &, F90geom &);
  void mpas_geo_info_f90(const F90geom &, int &, int &, int &, int &, int &,
                         int &, int &, int &);
  void mpas_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void mpas_model_setup_f90(const eckit::Configuration * const *,
                            const F90geom &, F90model &);
  void mpas_model_delete_f90(F90model &);

  void mpas_model_prepare_integration_f90(const F90model &, const F90flds &);
  void mpas_model_prepare_integration_tl_f90(const F90model &, const F90flds &);
  void mpas_model_prepare_integration_ad_f90(const F90model &, const F90flds &);

  void mpas_model_propagate_f90(const F90model &, const F90flds &);
  void mpas_model_prop_traj_f90(const F90model &, const F90flds &, F90traj &);
  void mpas_model_propagate_tl_f90(const F90model &, const F90flds &,
                                   const F90traj &);
  void mpas_model_propagate_ad_f90(const F90model &, const F90flds &,
                                   const F90traj &);

  void mpas_model_wipe_traj_f90(F90traj &);
  void mpas_traj_minmaxrms_f90(const F90traj &, double &);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void mpas_field_create_f90(F90flds &, const F90geom &,
                             const eckit::Configuration * const *);
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
  void mpas_field_diff_incr_f90(const F90flds &, const F90flds &,
                                const F90flds &);

  void mpas_field_change_resol_f90(const F90flds &, const F90flds &);

  void mpas_field_read_file_f90(const F90flds &,
                                const eckit::Configuration * const *,
                                util::DateTime * const *);
  void mpas_field_analytic_init_f90(const F90flds &, const F90geom &,
                                        const eckit::Configuration * const *,
                                        util::DateTime * const *);
  void mpas_field_write_file_f90(const F90flds &,
                                 const eckit::Configuration * const *,
                                 const util::DateTime * const *);

  void mpas_field_getvalues_notraj_f90(const F90flds &, const F90locs &,
                             const eckit::Configuration * const *,
                             const F90goms &);
  void mpas_field_getvalues_f90(const F90flds &, const F90locs &,
                             const eckit::Configuration * const *,
                             const F90goms &, const F90ootrj &);
  void mpas_field_getvalues_tl_f90(const F90flds &, const F90locs &,
                             const eckit::Configuration * const *,
                             const F90goms &, const F90ootrj &);
  void mpas_field_getvalues_ad_f90(const F90flds &, const F90locs &,
                             const eckit::Configuration * const *,
                             const F90goms &, const F90ootrj &);
  void mpas_field_ug_coord_f90(const F90flds &, const int &, const int &);
  void mpas_field_field_to_ug_f90(const F90flds &, const int &, const int &);
  void mpas_field_field_from_ug_f90(const F90flds &, const int &);

  void mpas_field_gpnorm_f90(const F90flds &, const int &, double &);
  void mpas_field_sizes_f90(const F90flds &, int &, int &);
  void mpas_field_rms_f90(const F90flds &, double &);

  void mpas_field_dirac_f90(const F90flds &,
                            const eckit::Configuration * const *);

  void mpas_getvaltraj_setup_f90(const F90ootrj &);
  void mpas_getvaltraj_delete_f90(const F90ootrj &);

// -----------------------------------------------------------------------------
//  Background error
// -----------------------------------------------------------------------------
  void mpas_b_setup_f90(F90bmat &, const eckit::Configuration * const *,
                        const F90geom &);
  void mpas_b_delete_f90(F90bmat &);

  void mpas_b_linearize_f90(const F90bmat &,
                            const eckit::Configuration * const *);

  void mpas_b_mult_f90(const F90bmat &, const F90flds &, const F90flds &);
  void mpas_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);

  void mpas_b_randomize_f90(const F90bmat &, const F90flds &);

// -----------------------------------------------------------------------------
//  Variable Change for Background error matrix
// -----------------------------------------------------------------------------
  void mpas_varchange_setup_f90(const F90vc &, const F90flds &,
                                   const F90flds &, const F90geom &,
                                   const eckit::Configuration * const *);
  void mpas_varchange_delete_f90(F90vc &);
  void mpas_varchange_multiply_f90(const F90vc &, const F90flds &,
                                      const F90flds &);
  void mpas_varchange_multiplyadjoint_f90(const F90vc &, const F90flds &,
                                      const F90flds &);
  void mpas_varchange_multiplyinverse_f90(const F90vc &, const F90flds &,
                                      const F90flds &);
  void mpas_varchange_multiplyinverseadjoint_f90(const F90vc &,
                                              const F90flds &, const F90flds &);
// -----------------------------------------------------------------------------
//  Localization matrix
// -----------------------------------------------------------------------------
  void mpas_localization_setup_f90(F90lclz &,
                                 const eckit::Configuration * const *,
                                 const F90geom &);
  void mpas_localization_delete_f90(F90lclz &);
  void mpas_localization_mult_f90(const F90lclz &, const F90flds &);

// -----------------------------------------------------------------------------
//  Run
// -----------------------------------------------------------------------------
  void mpas_init_f90(const eckit::Configuration * const *);
  void mpi_finalize_f90();

}
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_FORTRAN_H_
