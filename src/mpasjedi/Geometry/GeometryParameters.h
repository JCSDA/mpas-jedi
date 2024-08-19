/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

/// Configuration options recognized by mpas_geom_mod
class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  /// namelist filename to be used in mpas_init
  oops::RequiredParameter<std::string> nml_file{"nml_file", this};

  /// streams filename to be used in mpas_init
  oops::RequiredParameter<std::string> streams_file{"streams_file", this};

  /// option to deallocate not-used fields for reducing memory usage
  /// can not be true for forecast and hofx (4D) applications
  oops::Parameter<bool> deallocate_non_da_fields{ "deallocate non-da fields", false, this};

  /// yaml filename that contains the list of variables to be kept
  /// when "deallocate non-da fields" is true
  oops::Parameter<std::string> kept_fields_file{ "kept fields file", "keptvars.yaml", this};

  /// yaml filename that contains configurations of templated field names
  oops::Parameter<std::string> template_fields_file{ "template fields file", "geovars.yaml", this};

  /// vertical coordinate for BUMP to be used in the parameter estimate application
  /// not needed for hofx and variational applications
  /// available options are modellevel, height, avgheight, and scaleheight
  oops::Parameter<std::string> bump_vunit{ "bump vunit", "modellevel", this};

  /// LocalEnsembleDA iterator dimensionality (2 or 3)
  oops::Parameter<int> iterator_dimension{ "iterator dimension", 2, this};

  /// gmsh
  oops::Parameter<bool> gmsh_save{ "gmsh_save", false, this};
  oops::Parameter<std::string> gmsh_filename{ "gmsh_filename", "out.msh", this};

  /// SACA-related parameters
  /// NOTE:: "l_build_madwrf" and "l_build_gsdcloud" are mutually exclusive
  ///        "l_saturate_qv" and "l_conserve_thetaV" are mutually exclusive
  /// whether to use MADWRF's cloud building algorithm
  oops::Parameter<bool> l_build_madwrf{ "l_build_madwrf", true, this};
  /// whether to use GSD's cloud building algorithm
  oops::Parameter<bool> l_build_gsdcloud{ "l_build_gsdcloud", false, this};
  /// whether to saturate water vapor mixing ratio (Qv)
  oops::Parameter<bool> l_saturate_qv{ "l_saturate_qv", false, this};
  /// whether to conserve Theta V
  oops::Parameter<bool> l_conserve_thetaV{ "l_conserve_thetaV", true, this};
  /// cloud fraction value for cloud insertion
  oops::Parameter<float> cldfra_def{ "cldfra_def", 0.98f, this};
  /// cloud fraction threshold to determine presence of model cloud
  /// only applicable to l_build_madwrf=true
  oops::Parameter<float> cldfra_thresh{ "cldfra_thresh", 0.0f, this};
  /// threshold to distinguish clear/cloudy for "observed" cldmask (or cloud fraction)
  oops::Parameter<float> cldmask_thresh{ "cldmask_thresh", 0.0f, this};
  /// above ground height [meters] to limit the cloud building
  /// only works with l_build_gsdcloud=true
  oops::Parameter<float> cld_bld_hgt{ "cld_bld_hgt", 1200.0f, this};
};
// -------------------------------------------------------------------------------------------------

}  // namespace mpas
