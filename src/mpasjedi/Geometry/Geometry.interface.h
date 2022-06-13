/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include "mpasjedi/Fortran.h"

// Forward declarations
namespace atlas {
  namespace field {
    class FieldSetImpl;
  }
  namespace functionspace {
    class FunctionSpaceImpl;
  }
}

namespace eckit {
  class Configuration;
  namespace mpi {
    class Comm;
  }
}

namespace oops {
  class Variables;
}

namespace mpas {
extern "C" {
  void mpas_geo_setup_f90(F90geom &,
                          const eckit::Configuration &,
                          const eckit::mpi::Comm *);
  void mpas_geo_set_lonlat_f90(const F90geom &, atlas::field::FieldSetImpl *,
                                     const bool &);
  void mpas_geo_set_functionspace_pointer_f90(const F90geom &,
                                                    atlas::functionspace::FunctionSpaceImpl *,
                                                    atlas::functionspace::FunctionSpaceImpl *);
  void mpas_geo_fill_extra_fields_f90(const F90geom &,
                                      atlas::field::FieldSetImpl *);
  void mpas_geo_clone_f90(F90geom &, const F90geom &);
  void mpas_geo_is_equal_f90(bool &, const F90geom &, const F90geom &);
  void mpas_geo_vars_nlevels_f90(const F90geom &, const oops::Variables &,
                                 const std::size_t &, std::size_t &);
  void mpas_geo_info_f90(const F90geom &, int &, int &, int &, int &, int &,
                         int &, int &, int &);
  void mpas_geo_delete_f90(F90geom &);
}
// -----------------------------------------------------------------------------

}  // namespace mpas
