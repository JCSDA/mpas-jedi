/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

#include "mpasjedi/Geometry/Geometry.interface.h"
#include "mpasjedi/Geometry/GeometryParameters.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Geometry handles geometry for MPAS model.

class Geometry : public util::Printable,
                      private util::ObjectCounter<Geometry> {
 public:
  typedef GeometryParameters Parameters_;
  static const std::string classname() {return "mpas::Geometry";}

  Geometry(const Parameters_ &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FunctionSpace * atlasFunctionSpaceIncludingHalo() const {
      return atlasFunctionSpaceIncludingHalo_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  void latlon(std::vector<real_type> &, std::vector<real_type> &, const bool) const;

  bool isEqual(const Geometry &) const;

  std::vector<size_t> variableSizes(const oops::Variables &) const;

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpaceIncludingHalo_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
