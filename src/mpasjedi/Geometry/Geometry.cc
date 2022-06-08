/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
Geometry::Geometry(const GeometryParameters & params,
                           const eckit::mpi::Comm & comm) : comm_(comm) {
  oops::Log::trace() << "========= Geometry::Geometry step 1 =========="
                     << std::endl;
  mpas_geo_setup_f90(keyGeom_, params.toConfiguration(), &comm);

  // Set ATLAS lon/lat field
  atlasFieldSet_.reset(new atlas::FieldSet());
  const bool include_halo = true;  // include halo so we can set up both function spaces
  mpas_geo_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get(), include_halo);

  // Create ATLAS function space
  const atlas::Field atlasField = atlasFieldSet_->field("lonlat");
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));
  ASSERT(include_halo);
  // Create ATLAS function space with halo
  const atlas::Field atlasFieldInclHalo = atlasFieldSet_->field("lonlat_including_halo");
  atlasFunctionSpaceIncludingHalo_.reset(new atlas::functionspace::PointCloud(atlasFieldInclHalo));

  // Set ATLAS function space pointer in Fortran
  mpas_geo_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get(),
                                               atlasFunctionSpaceIncludingHalo_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  mpas_geo_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());

  oops::Log::trace() << "========= Geometry::Geometry step 2 =========="
                     << std::endl;
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  oops::Log::trace() << "========= Geometry mpas_geo_clone_f90   =========="
                     << std::endl;
  mpas_geo_clone_f90(keyGeom_, other.keyGeom_);

  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                                       other.atlasFunctionSpace_->lonlat()));
  atlasFunctionSpaceIncludingHalo_.reset(new atlas::functionspace::PointCloud(
                            other.atlasFunctionSpaceIncludingHalo_->lonlat()));
  mpas_geo_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get(),
                                               atlasFunctionSpaceIncludingHalo_->get());
  atlasFieldSet_.reset(new atlas::FieldSet());

  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
  mpas_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
bool Geometry::isEqual(const Geometry & other) const {
  bool isEqual;

  mpas_geo_is_equal_f90(isEqual, keyGeom_, other.keyGeom_);

  return isEqual;
}
// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  // vector of level counts
  std::vector<size_t> varSizes(vars.size());

  mpas_geo_vars_nlevels_f90(keyGeom_, vars, vars.size(), varSizes[0]);

  return varSizes;
}
// -----------------------------------------------------------------------------
void Geometry::latlon(std::vector<real_type> & lats, std::vector<real_type> & lons,
                      const bool halo) const {
  const atlas::functionspace::PointCloud * fspace;
  if (halo) {
    fspace = atlasFunctionSpaceIncludingHalo_.get();
  } else {
    fspace = atlasFunctionSpace_.get();
  }
  const auto lonlat = atlas::array::make_view<real_type, 2>(fspace->lonlat());
  const size_t npts = fspace->size();
  lats.resize(npts);
  lons.resize(npts);
  for (size_t jj = 0; jj < npts; ++jj) {
    lats[jj] = lonlat(jj, 1);
    lons[jj] = lonlat(jj, 0);
    if (lons[jj] < 0.0) lons[jj] += 360.0;
  }
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  int nCellsGlobal;
  int nCells;
  int nCellsSolve;
  int nEdgesGlobal;
  int nEdges;
  int nEdgesSolve;
  int nVertLevels;
  int nVertLevelsP1;
  mpas_geo_info_f90(keyGeom_, nCellsGlobal, nCells, nCellsSolve, \
                              nEdgesGlobal, nEdges, nEdgesSolve, \
                              nVertLevels, nVertLevelsP1);

  os << ", nCellsGlobal = " << nCellsGlobal \
     << ", nCells = " << nCells \
     << ", nCellsSolve = " << nCellsSolve \
     << ", nEdgesGlobal = " << nEdgesGlobal \
     << ", nEdges = " << nEdges \
     << ", nEdgesSolve = " << nEdgesSolve \
     << ", nVertLevels = " <<nVertLevels \
     << ", nVertLevelsP1 = " <<nVertLevelsP1
     << ", communicator = " << this->getComm().name();
}
// -----------------------------------------------------------------------------
}  // namespace mpas
