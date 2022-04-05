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

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_set_atlas_lonlat(c_key_self, c_afieldset, c_include_halo) bind(c,name='mpas_geo_set_atlas_lonlat_f90')
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

call geo_set_atlas_lonlat(self, afieldset, c_include_halo)

end subroutine c_mpas_geo_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_set_atlas_functionspace_pointer(c_key_self,c_afunctionspace,c_afunctionspace_incl_halo) &
 & bind(c,name='mpas_geo_set_atlas_functionspace_pointer_f90')
!!use atlas_module, only: atlas_functionspace_pointcloud
use atlas_module, only: atlas_fieldset, atlas_functionspace
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afunctionspace
type(c_ptr), intent(in), value :: c_afunctionspace_incl_halo
type(mpas_geom),pointer :: self

call mpas_geom_registry%get(c_key_self, self)
self%afunctionspace = atlas_functionspace(c_afunctionspace)
self%afunctionspace_incl_halo = atlas_functionspace(c_afunctionspace_incl_halo)

end subroutine c_mpas_geo_set_atlas_functionspace_pointer

! --------------------------------------------------------------------------------------------------

subroutine c_mpas_geo_fill_atlas_fieldset(c_key_self, c_afieldset) &
 & bind(c,name='mpas_geo_fill_atlas_fieldset_f90')
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

call geo_fill_atlas_fieldset(self, afieldset)

end subroutine c_mpas_geo_fill_atlas_fieldset

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
