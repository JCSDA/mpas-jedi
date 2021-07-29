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
  void mpas_geo_vars_nlevels_f90(const F90geom &, const oops::Variables &,
                                 const std::size_t &, std::size_t &);
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

}
// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_FORTRAN_H_
