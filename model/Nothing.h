/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPAS_MODEL_NOTHING_H_
#define MPAS_MODEL_NOTHING_H_

#include <ostream>

#include "oops/util/Printable.h"

namespace mpas {

// -----------------------------------------------------------------------------

class Nothing : public util::Printable {
 public:
  Nothing() {}
  ~Nothing() {}

 private:
  void print(std::ostream &) const {}
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPAS_MODEL_NOTHING_H_
