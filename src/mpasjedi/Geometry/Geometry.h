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
  const atlas::FunctionSpace & functionSpace() const {return functionSpaceIncludingHalo_;}
  atlas::FunctionSpace & functionSpace() {return functionSpaceIncludingHalo_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}
  void latlon(std::vector<real_type> &, std::vector<real_type> &, const bool) const;

  bool isEqual(const Geometry &) const;
  // The bool levelsAreTopDown=true is a fixed parameter here (cannot be false);
  // it is a required parameter to be passed to oops::Geometry -> getvalue -> ufo
  // - MPAS-Model levels are bottom-to-up, differs from CRTM: Top-to-Bottom
  // - In mpas-jedi, flip happens in toFieldSet/fromFieldSet(atlas)/geo_fill_extra_fields(saber)
  const bool levelsAreTopDown() const {return true;}

  std::vector<size_t> variableSizes(const oops::Variables &) const;

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace functionSpace_;
  atlas::FunctionSpace functionSpaceIncludingHalo_;
  atlas::FieldSet extraFields_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
