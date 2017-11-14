/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "model/Variables.h"
#include "model/Fortran.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  int nv;
  int nl;
  mpas_var_info_f90(keyVar_, nv, nl);
  os << nv;
  if (nl == 1) os << " with LBC";
  ASSERT(nl == 0 || nl == 1);
}

// -----------------------------------------------------------------------------

}  // namespace mpas
