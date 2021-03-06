/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPASJEDI_GEOMETRYMPAS_H_
#define MPASJEDI_GEOMETRYMPAS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/Fortran.h"

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

  bool isEqual(const GeometryMPAS &) const;

  std::vector<size_t> variableSizes(const oops::Variables &) const;

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

#endif  // MPASJEDI_GEOMETRYMPAS_H_
