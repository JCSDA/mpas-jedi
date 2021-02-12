/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPASJEDI_FORTRAN_H_
#define MPASJEDI_FORTRAN_H_

#include <atlas/field.h>
#include <atlas/functionspace.h>

#include "oops/base/Variables.h"

// Forward declarations
namespace eckit {
  class Configuration;
  namespace mpi {
    class Comm;
  }
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
// LinearVariableChange key
typedef int F90lvcc2a;
// VarChange key
typedef int F90vc;
// GetValues key
typedef int F90getvalues;
// LinearGetValues key
typedef int F90lineargetvalues;

/// Interface to Fortran MPAS model
/*!
 * The core of the MPAS model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Run
// -----------------------------------------------------------------------------
  void mpas_run_init_f90(const eckit::Configuration &);
  void mpas_run_final_f90();

// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mpas_geo_setup_f90(F90geom &,
                          const eckit::Configuration &,
                          const eckit::mpi::Comm *);
  void mpas_geo_set_atlas_lonlat_f90(const F90geom &,
                                     atlas::field::FieldSetImpl *);
  void mpas_geo_set_atlas_functionspace_pointer_f90(const F90geom &,
                    atlas::functionspace::FunctionSpaceImpl *);
  void mpas_geo_fill_atlas_fieldset_f90(const F90geom &,
                                        atlas::field::FieldSetImpl *);
  void mpas_geo_clone_f90(F90geom &, const F90geom &);
  void mpas_geo_is_equal_f90(bool &, const F90geom &, const F90geom &);
  void mpas_geo_info_f90(const F90geom &, int &, int &, int &, int &, int &,
                         int &, int &, int &);
  void mpas_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
//  ErrorCovariance (Background error)
// -----------------------------------------------------------------------------
  void mpas_b_setup_f90(F90bmat &, const eckit::Configuration &,
                        const F90geom &);
  void mpas_b_delete_f90(F90bmat &);

  void mpas_b_linearize_f90(const F90bmat &,
                            const eckit::Configuration &);

  void mpas_b_mult_f90(const F90bmat &, const F90inc &, const F90inc &);
  void mpas_b_invmult_f90(const F90bmat &, const F90inc &, const F90inc &);

  void mpas_b_randomize_f90(const F90bmat &, const F90inc &);

// -----------------------------------------------------------------------------
//  VarChange (Variable Change for Background error matrix)
// -----------------------------------------------------------------------------
  void mpas_varchange_setup_f90(const F90vc &, const F90state &,
                                   const F90state &, const F90geom &,
                                   const eckit::Configuration &);
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
//  LinearVariableChange (Linear variable Change for Background error matrix)
// -----------------------------------------------------------------------------
  void mpas_linvarcha_c2a_setup_f90(const F90lvcc2a &, const F90state &,
                                   const F90state &, const F90geom &,
                                   const eckit::Configuration &, const oops::Variables &);
  void mpas_linvarcha_c2a_delete_f90(F90lvcc2a &);
  void mpas_linvarcha_c2a_multiply_f90(const F90lvcc2a &, const F90inc &,
                                      const F90inc &);
  void mpas_linvarcha_c2a_multiplyadjoint_f90(const F90lvcc2a &, const F90inc &,
                                      const F90inc &);
  void mpas_linvarcha_c2a_multiplyinverse_f90(const F90lvcc2a &, const F90inc &,
                                      const F90inc &);
  void mpas_linvarcha_c2a_multiplyinverseadjoint_f90(const F90lvcc2a &,
                                              const F90inc &, const F90inc &);
// -----------------------------------------------------------------------------
//  LocalizationMatrix
// -----------------------------------------------------------------------------
  void mpas_localization_setup_f90(F90lclz &,
                                 const eckit::Configuration &,
                                 const F90geom &);
  void mpas_localization_delete_f90(F90lclz &);
  void mpas_localization_randomize_f90(const F90lclz &, const F90inc &);
  void mpas_localization_mult_f90(const F90lclz &, const F90inc &);

}
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_FORTRAN_H_
