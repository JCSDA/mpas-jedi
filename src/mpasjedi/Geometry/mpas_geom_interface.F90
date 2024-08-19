! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------
subroutine c_mpas_geo_setup(c_key_self, c_conf, c_comm) bind(c,name='mpas_geo_setup_f90')
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_comm

type(mpas_geom), pointer :: self
type(fckit_configuration) :: f_conf

call mpas_geom_registry%init()
call mpas_geom_registry%add(c_key_self)
call mpas_geom_registry%get(c_key_self, self)

f_conf = fckit_configuration(c_conf)
call geo_setup(self, f_conf, fckit_mpi_comm(c_comm))

end subroutine c_mpas_geo_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_clone(c_key_self, c_key_other) bind(c,name='mpas_geo_clone_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(inout)    :: c_key_self
integer(c_int), intent(in) :: c_key_other

type(mpas_geom), pointer :: self, other

call mpas_geom_registry%add(c_key_self)
call mpas_geom_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_other, other)

call geo_clone(self, other)

end subroutine c_mpas_geo_clone

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_delete(c_key_self) bind(c,name='mpas_geo_delete_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(inout) :: c_key_self     
type(mpas_geom), pointer :: self

call mpas_geom_registry%get(c_key_self, self)
call geo_delete(self)
call mpas_geom_registry%remove(c_key_self)

end subroutine c_mpas_geo_delete

! ------------------------------------------------------------------------------

!> C++ interface to get dimension of the GeometryIterator
subroutine mpas_geo_iterator_dimension(c_key_self, itd) bind(c, name='mpas_geo_iterator_dimension_f90')
  use iso_c_binding
  use mpas_geom_mod

  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(out) :: itd ! iterator dimension

  type(mpas_geom), pointer :: self
  call mpas_geom_registry%get(c_key_self, self)

  itd = self%iterator_dimension
end subroutine mpas_geo_iterator_dimension

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_set_lonlat(c_key_self, c_afieldset, c_include_halo) bind(c,name='mpas_geo_set_lonlat_f90')
use atlas_module, only: atlas_fieldset
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in), value :: c_afieldset
logical, intent(in) :: c_include_halo
type(mpas_geom), pointer :: self
type(atlas_fieldset) :: afieldset

call mpas_geom_registry%get(c_key_self, self)
afieldset = atlas_fieldset(c_afieldset)

call geo_set_lonlat(self, afieldset, c_include_halo)

end subroutine c_mpas_geo_set_lonlat

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_set_functionspace_pointer(c_key_self,c_afunctionspace) &
 & bind(c,name='mpas_geo_set_functionspace_pointer_f90')
use atlas_module, only: atlas_fieldset, atlas_functionspace
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afunctionspace
type(mpas_geom),pointer :: self

call mpas_geom_registry%get(c_key_self, self)
self%afunctionspace = atlas_functionspace(c_afunctionspace)

end subroutine c_mpas_geo_set_functionspace_pointer

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_fill_geometry_fields(c_key_self, c_afieldset) &
 & bind(c,name='mpas_geo_fill_geometry_fields_f90')
use atlas_module, only: atlas_fieldset
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_afieldset
type(mpas_geom), pointer :: self
type(atlas_fieldset) :: afieldset

call mpas_geom_registry%get(c_key_self, self)
afieldset = atlas_fieldset(c_afieldset)

call geo_fill_geometry_fields(self, afieldset)

end subroutine c_mpas_geo_fill_geometry_fields

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_is_equal(c_is_equal, c_key_self, c_key_other) &
 & bind(c,name='mpas_geo_is_equal_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
logical(c_bool),    intent(inout) :: c_is_equal
integer(c_int),     intent(in)    :: c_key_self
integer(c_int),     intent(in)    :: c_key_other
type(mpas_geom), pointer :: self
type(mpas_geom), pointer :: other

call mpas_geom_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_other, other)

call geo_is_equal(c_is_equal, self, other)

end subroutine c_mpas_geo_is_equal

! ------------------------------------------------------------------------------

subroutine mpas_geo_vars_nlevels_c(c_key_self, c_vars, c_nvars, c_nlevels) &
      bind(c,name='mpas_geo_vars_nlevels_f90')
use iso_c_binding
use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
use oops_variables_mod, only: oops_variables
implicit none
integer(c_int), intent(in) :: c_key_self !< Geometry
type(c_ptr), value, intent(in) :: c_vars !< Variables
integer(c_size_t), intent(in) :: c_nvars !< size of Variables
integer(c_size_t), intent(inout) :: c_nlevels(c_nvars) !< level counts

type(mpas_geom), pointer :: self
type(oops_variables) :: vars

call mpas_geom_registry%get(c_key_self, self)
vars = oops_variables(c_vars)
call self%nlevels(vars, c_nvars, c_nlevels)

end subroutine mpas_geo_vars_nlevels_c

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_info(c_key_self, c_nCellsGlobal, c_nCells, c_nCellsSolve, &
                                       c_nEdgesGlobal, c_nEdges, c_nEdgesSolve, &
                                       c_nVertLevels, c_nVertLevelsP1) &
                                       bind(c,name='mpas_geo_info_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: &
   c_nCellsGlobal, c_nCells, c_nCellsSolve, &
   c_nEdgesGlobal, c_nEdges, c_nEdgesSolve, &
   c_nVertLevels, c_nVertLevelsP1

type(mpas_geom), pointer :: self

call mpas_geom_registry%get(c_key_self, self)
call geo_info(self, c_nCellsGlobal, c_nCells, c_nCellsSolve, &
                    c_nEdgesGlobal, c_nEdges, c_nEdgesSolve, &
                    c_nVertLevels, c_nVertLevelsP1)

end subroutine c_mpas_geo_info

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_vert_coord(c_key_self, c_nVertLevels, c_CoordNameLen, &
                                 c_CoordName, c_VertCoord) &
                                 bind(c,name='mpas_geo_vert_coord_f90')
use iso_c_binding
use mpas_geom_mod
use mpas_kind_types
use mpas_kinds, only : c_real_type
implicit none
integer(c_int), intent(in)       :: c_key_self
integer(c_int), intent(in)       :: c_nVertLevels
integer(c_int), intent(in)       :: c_CoordNameLen
character(c_char), intent(in)    :: c_CoordName(c_CoordNameLen)
real(c_real_type), intent(out)   :: c_VertCoord(c_nVertLevels)

type(mpas_geom), pointer :: self
character(len=c_CoordNameLen) :: coordname

! convert c string to fortran string
coordname = transfer(c_CoordName,coordname)

call mpas_geom_registry%get(c_key_self, self)

call geo_vert_coord(self, c_nVertLevels, coordname, c_VertCoord)

end subroutine c_mpas_geo_vert_coord

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_get_num_nodes_and_elements(c_key_self, c_num_nodes, c_num_tris) &
    bind(c, name='mpas_geo_get_num_nodes_and_elements_f90')
use iso_c_binding
use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
implicit none
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: c_num_nodes
  integer(c_int), intent(out) :: c_num_tris

  integer :: num_nodes
  integer :: num_tris

  type(mpas_geom), pointer :: self
  call mpas_geom_registry%get(c_key_self, self)
  call self%get_num_nodes_and_elements(num_nodes, num_tris)

  c_num_nodes = num_nodes
  c_num_tris = num_tris

end subroutine c_mpas_geo_get_num_nodes_and_elements

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_get_coords_and_connectivities(c_key_self, &
    c_num_nodes, c_lons, c_lats, c_ghosts, c_global_indices, c_remote_indices, c_partition, &
    c_num_tri_boundary_nodes, c_raw_tri_boundary_nodes) &
    bind(c, name='mpas_geo_get_coords_and_connectivities_f90')
use iso_c_binding
use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
use kinds, only: kind_real !oops
implicit none
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent( in) :: c_num_nodes
  real(c_double), intent(out) :: c_lons(c_num_nodes)
  real(c_double), intent(out) :: c_lats(c_num_nodes)
  integer(c_int), intent(out) :: c_ghosts(c_num_nodes)
  integer(c_int), intent(out) :: c_global_indices(c_num_nodes)
  integer(c_int), intent(out) :: c_remote_indices(c_num_nodes)
  integer(c_int), intent(out) :: c_partition(c_num_nodes)
  integer(c_int), intent( in) :: c_num_tri_boundary_nodes
  integer(c_int), intent(out) :: c_raw_tri_boundary_nodes(c_num_tri_boundary_nodes)

  integer :: num_nodes, num_tri_boundary_nodes
  type(mpas_geom), pointer :: self

  real(kind_real) :: lons(c_num_nodes)
  real(kind_real) :: lats(c_num_nodes)
  integer :: ghosts(c_num_nodes)
  integer :: global_indices(c_num_nodes)
  integer :: remote_indices(c_num_nodes)
  integer :: partition(c_num_nodes)
  integer :: raw_tri_boundary_nodes(c_num_tri_boundary_nodes)

  num_nodes = c_num_nodes
  num_tri_boundary_nodes = c_num_tri_boundary_nodes

  call mpas_geom_registry%get(c_key_self, self)
  call self%get_coords_and_connectivities(num_nodes, num_tri_boundary_nodes, &
    lons, lats, ghosts, global_indices, remote_indices, partition, raw_tri_boundary_nodes)

  c_lons = lons
  c_lats = lats
  c_ghosts = ghosts
  c_global_indices = global_indices
  c_remote_indices = remote_indices
  c_partition = partition
  c_raw_tri_boundary_nodes = raw_tri_boundary_nodes

end subroutine c_mpas_geo_get_coords_and_connectivities

! ------------------------------------------------------------------------------
