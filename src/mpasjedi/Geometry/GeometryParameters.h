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

  /// interpolation method to be used to interpolate data between different geometries
  /// available options are bump and unstructured
  /// bump is more memory intensive than unstructured
  oops::Parameter<std::string> interpolation_type{ "interpolation type", "unstructured", this};

  /// option to deallocate not-used fields for reducing memory usage
  /// can not be true for forecast and hofx (4D) applications
  oops::Parameter<bool> deallocate_non_da_fields{ "deallocate non-da fields", false, this};

  /// yaml filename that contains configurations of templated field names
  oops::Parameter<std::string> template_fields_file{ "template fields file", "geovars.yaml", this};

  /// vertical coordinate for BUMP to be used in the parameter estimate application
  /// not needed for hofx and variational applications
  /// available options are modellevel, height, avgheight, and scaleheight
  oops::Parameter<std::string> bump_vunit{ "bump vunit", "modellevel", this};
};
// -------------------------------------------------------------------------------------------------

}  // namespace mpas
