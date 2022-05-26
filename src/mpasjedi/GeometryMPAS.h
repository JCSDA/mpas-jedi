/*
 * (C) Copyright 2017-2021 UCAR
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
#include "mpasjedi/GeometryMPASParameters.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// GeometryMPAS handles geometry for MPAS model.

class GeometryMPAS : public util::Printable,
                      private util::ObjectCounter<GeometryMPAS> {
 public:
  typedef GeometryMPASParameters Parameters_;
  static const std::string classname() {return "mpas::GeometryMPAS";}

  GeometryMPAS(const Parameters_ &, const eckit::mpi::Comm &);
  GeometryMPAS(const GeometryMPAS &);
  ~GeometryMPAS();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FunctionSpace * atlasFunctionSpaceIncludingHalo() const {
      return atlasFunctionSpaceIncludingHalo_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  void latlon(std::vector<real_type> &, std::vector<real_type> &, const bool) const;

  bool isEqual(const GeometryMPAS &) const;

  std::vector<size_t> variableSizes(const oops::Variables &) const;

 private:
  GeometryMPAS & operator=(const GeometryMPAS &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpaceIncludingHalo_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas

#endif  // MPASJEDI_GEOMETRYMPAS_H_
