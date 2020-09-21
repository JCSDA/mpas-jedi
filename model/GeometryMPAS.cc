/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"
#include "model/GeometryMPAS.h"
#include "Fortran.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
GeometryMPAS::GeometryMPAS(const eckit::Configuration & config,
                           const eckit::mpi::Comm & comm) : comm_(comm) {
  oops::Log::trace() << "========= GeometryMPAS::GeometryMPAS step 1 =========="
                     << std::endl;
  mpas_geo_setup_f90(keyGeom_, config, &comm);

  // Set ATLAS lon/lat field
  atlasFieldSet_.reset(new atlas::FieldSet());
  mpas_geo_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
  atlas::Field atlasField = atlasFieldSet_->field("lonlat");

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

  // Set ATLAS function space pointer in Fortran
  mpas_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
    atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  mpas_geo_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());

  oops::Log::trace() << "========= GeometryMPAS::GeometryMPAS step 2 =========="
                     << std::endl;
}
// -----------------------------------------------------------------------------
GeometryMPAS::GeometryMPAS(const GeometryMPAS & other) : comm_(other.comm_) {
  oops::Log::trace() << "========= GeometryMPAS mpas_geo_clone_f90   =========="
                     << std::endl;
  mpas_geo_clone_f90(keyGeom_, other.keyGeom_);
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                              other.atlasFunctionSpace_->lonlat()));
  mpas_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
    atlasFunctionSpace_->get());
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}
// -----------------------------------------------------------------------------
GeometryMPAS::~GeometryMPAS() {
  mpas_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
void GeometryMPAS::print(std::ostream & os) const {
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
