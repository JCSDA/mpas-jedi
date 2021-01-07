/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mpasjedi/Fortran.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace ufo {
  class Locations;
}

namespace mpas {

extern "C" {

  void mpas_lineargetvalues_create_f90(F90lineargetvalues &, const F90geom &, const ufo::Locations &);

  void mpas_lineargetvalues_delete_f90(F90lineargetvalues &);

  void mpas_lineargetvalues_set_trajectory_f90(const F90lineargetvalues &, const F90geom &,
                                               const F90state &, const util::DateTime &,
                                               const util::DateTime &, const ufo::Locations &,
                                               const F90goms &);

  void mpas_lineargetvalues_fill_geovals_tl_f90(const F90lineargetvalues &, const F90geom &,
                                                const F90inc &, const util::DateTime &,
                                                const util::DateTime &, const ufo::Locations &,
                                                const F90goms &);

  void mpas_lineargetvalues_fill_geovals_ad_f90(const F90lineargetvalues &, const F90geom &,
                                                const F90inc &, const util::DateTime &,
                                                const util::DateTime &, const ufo::Locations &,
                                                const F90goms &);

};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
