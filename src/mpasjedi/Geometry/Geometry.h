/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mpasjedi/GeometryIterator/GeometryIterator.h"

#include "mpasjedi/Geometry/Geometry.interface.h"

namespace mpas {

// -----------------------------------------------------------------------------
/// Geometry handles geometry for MPAS model.

class Geometry : public util::Printable,
                      private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "mpas::Geometry";}

  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  const atlas::FieldSet & fields() const {return fields_;}

  bool isEqual(const Geometry &) const;
  // The bool levelsAreTopDown=true is a fixed parameter here (cannot be false);
  // it is a required parameter to be passed to oops::Geometry -> getvalue -> ufo
  // - MPAS-Model levels are bottom-to-up, differs from CRTM: Top-to-Bottom
  // - In mpas-jedi, flip happens in toFieldSet/fromFieldSet(atlas)/geo_fill_extra_fields(saber)
  bool levelsAreTopDown() const {return true;}

  GeometryIterator begin() const;
  GeometryIterator end() const;
  int IteratorDimension() const;
  std::vector<int> nIterLevs(const oops::Variables &) const;
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  std::vector<real_type> verticalCoord(std::string &) const;
  int getDim(const std::string &) const;

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace functionSpace_;
  atlas::FieldSet fields_;
};
// -----------------------------------------------------------------------------

}  // namespace mpas
