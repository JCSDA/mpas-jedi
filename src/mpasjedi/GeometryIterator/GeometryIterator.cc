/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/GeometryIterator/GeometryIterator.h"
#include "mpasjedi/GeometryIterator/GeometryIterator.interface.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace mpas {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry & geom,
                                   const int & cellIndex,
                                   const int & levIndex)
  : geom_(geom) {
  mpas_geom_iter_setup_f90(keyIter_, geom_.toFortran(), cellIndex, levIndex);
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator & other)
  : geom_(other.geom_) {
  mpas_geom_iter_clone_f90(keyIter_, other.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  mpas_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  bool equals;
  mpas_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals && geom_.isEqual(other.geom_));
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  bool equals;
  mpas_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return !equals;
}

// -----------------------------------------------------------------------------

eckit::geometry::Point3 GeometryIterator::operator*() const {
  double lat, lon, height;
  mpas_geom_iter_current_f90(keyIter_, lon, lat, height);
  return eckit::geometry::Point3(lon, lat, height);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  mpas_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

int GeometryIterator::iteratorDimension() const {
  int dimension;
  mpas_geom_iter_dimension_f90(keyIter_, dimension);
  return dimension;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  double lat, lon, height;
  mpas_geom_iter_current_f90(keyIter_, lon, lat, height);
  os << "GeometryIterator, lon/lat/height: " << lon << " / " << lat
     << " / " << height << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace mpas
