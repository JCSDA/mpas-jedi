/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_GETVALUESTRAJMPAS_H_
#define MODEL_GETVALUESTRAJMPAS_H_

#include <ostream>

#include "oops/util/Printable.h"
#include "Fortran.h"

namespace mpas {

// -----------------------------------------------------------------------------

class GetValuesTrajMPAS : public util::Printable {
 public:
  GetValuesTrajMPAS();
  ~GetValuesTrajMPAS();

  int & toFortran() {return keyGetValuesTraj_;}
  const int & toFortran() const {return keyGetValuesTraj_;}

 private:
  void print(std::ostream &) const {}
  F90ootrj keyGetValuesTraj_;
};

// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_GETVALUESTRAJMPAS_H_
