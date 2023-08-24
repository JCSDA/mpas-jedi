/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Geometry/GeometryParameters.h"

// -----------------------------------------------------------------------------
namespace mpas {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & config,
                           const eckit::mpi::Comm & comm) : comm_(comm) {
  oops::Log::trace() << "========= Geometry::Geometry step 1 =========="
                     << std::endl;
  GeometryParameters params;
  params.deserialize(config);
  mpas_geo_setup_f90(keyGeom_, params.toConfiguration(), &comm);

  // Set ATLAS lon/lat field with halo
  atlas::FieldSet fs;
  const bool include_halo = true;
  mpas_geo_set_lonlat_f90(keyGeom_, fs.get(), include_halo);

  // Create function space
  const atlas::Field lonlatField = fs.field("lonlat");
  functionSpace_ = atlas::functionspace::PointCloud(lonlatField);

  // Create function space with halo
  const atlas::Field lonlatFieldInclHalo = fs.field("lonlat_including_halo");
  functionSpaceIncludingHalo_ = atlas::functionspace::PointCloud(lonlatFieldInclHalo);

  // Set ATLAS function space pointer in Fortran
  mpas_geo_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get(),
                                             functionSpaceIncludingHalo_.get());

  // Fill extra fieldset : for saber vunit
  extraFields_ =  atlas::FieldSet();
  mpas_geo_fill_extra_fields_f90(keyGeom_, extraFields_.get());

  oops::Log::trace() << "========= Geometry::Geometry step 2 =========="
                     << std::endl;
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  oops::Log::trace() << "========= Geometry mpas_geo_clone_f90   =========="
                     << std::endl;
  mpas_geo_clone_f90(keyGeom_, other.keyGeom_);
  functionSpace_ = atlas::functionspace::PointCloud(other.functionSpace_.lonlat());
  functionSpaceIncludingHalo_ = atlas::functionspace::PointCloud(
    other.functionSpaceIncludingHalo_.lonlat());
  mpas_geo_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get(),
                                         functionSpaceIncludingHalo_.get());
  extraFields_ = atlas::FieldSet();
  for (auto & field : other.extraFields_) {
    extraFields_->add(field);
  }
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
  mpas_geo_delete_f90(keyGeom_);
}
// -------------------------------------------------------------------------------------------------
GeometryIterator Geometry::begin() const {
  return GeometryIterator(*this, 1, 1);
}

// -------------------------------------------------------------------------------------------------
GeometryIterator Geometry::end() const {
  // return end of the geometry on this mpi tile
  // (returns index out of bounds for the iterator loops to work)
  return GeometryIterator(*this, -1, -1);
}
// -----------------------------------------------------------------------------
int Geometry::IteratorDimension() const {
  // return dimesnion of the iterator
  // if 2, iterator is over vertical columns
  // if 3, iterator is over 3D points
  int rv;
  mpas_geo_iterator_dimension_f90(keyGeom_, rv);
  return rv;
}
// -----------------------------------------------------------------------------
std::vector<int> Geometry::nIterLevs(const oops::Variables & vars) const {
  std::vector<int> n(vars.size());
  if (IteratorDimension() == 3) {
    // 1 level per variable
    std::fill(n.begin(), n.end(), 1);
  } else {
    // nVertLevels for all variables (including 2d)
    std::fill(n.begin(), n.end(), getDim("nVertLevels"));
  }
  return n;
}
// -------------------------------------------------------------------------------------------------
std::vector<double> Geometry::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  // TODO(JJG): get this to work for height and scale height
  std::vector<double> vc(getDim("nVertLevels"));
  if (vcUnits == "modellevel") {
    for (size_t i=0; i < vc.size(); ++i) {vc[i] = static_cast<double> (i+1);}
  } else {
    std::stringstream errorMsg;
    errorMsg << "Uknown vertical coordinate unit " << vcUnits << std::endl;
    ABORT(errorMsg.str());
  }
  return vc;
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
  const atlas::FunctionSpace * fspace;
  if (halo) {
    fspace = &functionSpaceIncludingHalo_;
  } else {
    fspace = &functionSpace_;
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
  mpas_geo_info_f90(keyGeom_, nCellsGlobal, nCells, nCellsSolve,
                              nEdgesGlobal, nEdges, nEdgesSolve,
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
const int Geometry::getDim(const std::string & dim) const {
  int ncGlobal;
  int nc;
  int ncSolve;
  int neGlobal;
  int ne;
  int neSolve;
  int nv;
  int nvP1;
  mpas_geo_info_f90(keyGeom_, ncGlobal, nc, ncSolve,
                              neGlobal, ne, neSolve,
                              nv, nvP1);
  if (dim == "nCellsGlobal") {
      return ncGlobal;
  } else if (dim == "nCells") {
      return nc;
  } else if (dim == "nCellsSolve") {
      return ncSolve;
  } else if (dim == "nEdgesGlobal") {
      return neGlobal;
  } else if (dim == "nEdges") {
      return ne;
  } else if (dim == "nEdgesSolve") {
      return neSolve;
  } else if (dim == "nVertLevels") {
      return nv;
  } else if (dim == "nVertLevelsP1") {
      return nvP1;
  } else {
      std::stringstream errorMsg;
      errorMsg << "Uknown dimension " << dim << std::endl;
      ABORT(errorMsg.str());
  }
  return -1;
}
// -----------------------------------------------------------------------------
}  // namespace mpas
