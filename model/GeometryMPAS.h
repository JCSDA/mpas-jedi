/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_GEOMETRYMPAS_H_
#define MPAS_MODEL_GEOMETRYMPAS_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mpas {

// -----------------------------------------------------------------------------
/// GeometryMPAS handles geometry for MPAS model.

class GeometryMPAS : public util::Printable,
                      private util::ObjectCounter<GeometryMPAS> {
 public:
  static const std::string classname() {return "mpas::GeometryMPAS";}

  explicit GeometryMPAS(const eckit::Configuration &);
  GeometryMPAS(const GeometryMPAS &);
  ~GeometryMPAS();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}

 private:
  GeometryMPAS & operator=(const GeometryMPAS &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPAS_MODEL_GEOMETRYMPAS_H_
