/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_GEOMETRYMPAS_H_
#define MODEL_GEOMETRYMPAS_H_

#include <memory>
#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

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

  GeometryMPAS(const eckit::Configuration &, const eckit::mpi::Comm &);
  GeometryMPAS(const GeometryMPAS &);
  ~GeometryMPAS();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}

  atlas::FunctionSpace * atlasFunctionSpace() const
    {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const
   {return atlasFieldSet_.get();}


 private:
  GeometryMPAS & operator=(const GeometryMPAS &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MODEL_GEOMETRYMPAS_H_
