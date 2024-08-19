/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/output/Gmsh.h"

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

  // setup the atlas functionspace
  {
    // get the number of nodes and cells owned by this PE
    int num_nodes;
    int num_tri_elements;
    mpas_geo_get_num_nodes_and_elements_f90(keyGeom_, num_nodes, num_tri_elements);

    std::vector<double> lons(num_nodes);
    std::vector<double> lats(num_nodes);
    std::vector<int> ghosts(num_nodes);
    std::vector<int> global_indices(num_nodes);
    std::vector<int> remote_indices(num_nodes);
    std::vector<int> partitions(num_nodes);

    const int num_tri_nodes = 3*num_tri_elements;
    std::vector<int> raw_tri_boundary_nodes(num_tri_nodes);
    mpas_geo_get_coords_and_connectivities_f90(keyGeom_,
      num_nodes, lons.data(), lats.data(), ghosts.data(),
      global_indices.data(), remote_indices.data(), partitions.data(),
      num_tri_nodes, raw_tri_boundary_nodes.data());

    // calculate per-PE global tri numbering offset
    std::vector<int> num_elements_per_rank(comm_.size());
    comm_.allGather(num_tri_elements, num_elements_per_rank.begin(), num_elements_per_rank.end());
    int global_element_index = 0;
    for (size_t i = 0; i < comm_.rank(); ++i) {
      global_element_index += num_elements_per_rank[i];
    }

    // convert some of the temporary arrays into a form atlas expects
    using atlas::gidx_t;
    using atlas::idx_t;

    std::vector<gidx_t> atlas_global_indices(num_nodes);
    std::transform(global_indices.begin(), global_indices.end(), atlas_global_indices.begin(),
      [](const int index) {return atlas::gidx_t{index};});

    const atlas::idx_t remote_index_base = 1;  // 1-based indexing from Fortran
    std::vector<idx_t> atlas_remote_indices(num_nodes);
    std::transform(remote_indices.begin(), remote_indices.end(), atlas_remote_indices.begin(),
      [](const int index) {return atlas::idx_t{index};});

    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes(num_tri_elements);
    std::vector<gidx_t> tri_global_indices(num_tri_elements);
    for (size_t tri = 0; tri < num_tri_elements; ++tri) {
      for (size_t i = 0; i < 3; ++i) {
        tri_boundary_nodes[tri][i] = raw_tri_boundary_nodes[3*tri + i];
      }
      tri_global_indices[tri] = global_element_index++;  // work for both global/regional MPAS mesh
    }
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes{};  // MPAS does not have quad.
    std::vector<gidx_t> quad_global_indices{};  // MPAS does not have quad.

    // build the mesh!
    eckit::LocalConfiguration config{};
    config.set("mpi_comm", comm_.name());
    const atlas::mesh::MeshBuilder mesh_builder{};
    atlas::Mesh mesh = mesh_builder(
      lons, lats, ghosts,
      atlas_global_indices, atlas_remote_indices, remote_index_base, partitions,
      tri_boundary_nodes, tri_global_indices,
      quad_boundary_nodes, quad_global_indices, config);
    atlas::mesh::actions::build_halo(mesh, 1);
    functionSpace_ = atlas::functionspace::NodeColumns(mesh, config);

    // optionally save output for viewing with gmsh
    if (params.gmsh_save.value()) {
      std::string filename = params.gmsh_filename.value();
      atlas::output::Gmsh gmsh(filename,
          atlas::util::Config("coordinates", "xyz")
          | atlas::util::Config("ghost", true));  // enables viewing halos per task
      gmsh.write(mesh);
    }
  }

  // Set ATLAS lon/lat field with halo
  atlas::FieldSet fs;
  const bool include_halo = true;
  mpas_geo_set_lonlat_f90(keyGeom_, fs.get(), include_halo);

  // Set ATLAS function space pointer in Fortran
  mpas_geo_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get());

  // Fill geometry fieldset : for saber vunit
  fields_ =  atlas::FieldSet();
  mpas_geo_fill_geometry_fields_f90(keyGeom_, fields_.get());

  oops::Log::trace() << "========= Geometry::Geometry step 2 =========="
                     << std::endl;
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  oops::Log::trace() << "========= Geometry mpas_geo_clone_f90   =========="
                     << std::endl;
  mpas_geo_clone_f90(keyGeom_, other.keyGeom_);
  functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
  mpas_geo_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get());
  fields_ = atlas::FieldSet();
  for (auto & field : other.fields_) {
    fields_->add(field);
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
std::vector<real_type> Geometry::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  std::vector<real_type> vc(getDim("nVertLevels"));
  // vertical levels
  int vert_levels = vc.size();
  // length of vcUnits
  int len_vcunits = vcUnits.size();
  // vector of vcUnits
  std::vector<char> vec_vcUnits(vcUnits.begin(), vcUnits.end());
  vec_vcUnits.push_back('\0');  // Ensure null-termination for C compatibility.
  // vertical coordinates
  mpas_geo_vert_coord_f90(keyGeom_, vert_levels, len_vcunits, &vec_vcUnits[0], vc[0]);

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
int Geometry::getDim(const std::string & dim) const {
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
