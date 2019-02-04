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
// State key type
typedef int F90state;
// Increment key type
typedef int F90inc;
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

// -----------------------------------------------------------------------------
//  GetValuesTraj
// -----------------------------------------------------------------------------
  void mpas_getvaltraj_setup_f90(const F90ootrj &);
  void mpas_getvaltraj_delete_f90(const F90ootrj &);

// -----------------------------------------------------------------------------
//  ErrorCovariance (Background error)
// -----------------------------------------------------------------------------
  void mpas_b_setup_f90(F90bmat &, const eckit::Configuration * const *,
                        const F90geom &);
  void mpas_b_delete_f90(F90bmat &);

  void mpas_b_linearize_f90(const F90bmat &,
                            const eckit::Configuration * const *);

  void mpas_b_mult_f90(const F90bmat &, const F90inc &, const F90inc &);
  void mpas_b_invmult_f90(const F90bmat &, const F90inc &, const F90inc &);

  void mpas_b_randomize_f90(const F90bmat &, const F90inc &);

// -----------------------------------------------------------------------------
//  VarChange (Variable Change for Background error matrix)
// -----------------------------------------------------------------------------
  void mpas_varchange_setup_f90(const F90vc &, const F90state &,
                                   const F90state &, const F90geom &,
                                   const eckit::Configuration * const *);
  void mpas_varchange_delete_f90(F90vc &);
  void mpas_varchange_multiply_f90(const F90vc &, const F90inc &,
                                      const F90inc &);
  void mpas_varchange_multiplyadjoint_f90(const F90vc &, const F90inc &,
                                      const F90inc &);
  void mpas_varchange_multiplyinverse_f90(const F90vc &, const F90inc &,
                                      const F90inc &);
  void mpas_varchange_multiplyinverseadjoint_f90(const F90vc &,
                                              const F90inc &, const F90inc &);
// -----------------------------------------------------------------------------
//  LocalizationMatrix
// -----------------------------------------------------------------------------
  void mpas_localization_setup_f90(F90lclz &,
                                 const eckit::Configuration * const *,
                                 const F90geom &);
  void mpas_localization_delete_f90(F90lclz &);
  void mpas_localization_mult_f90(const F90lclz &, const F90inc &);

}
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MODEL_FORTRAN_H_
