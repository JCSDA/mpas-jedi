! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------
subroutine c_mpas_geo_setup(c_key_self, c_conf) bind(c,name='mpas_geo_setup_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in) :: c_conf

type(mpas_geom), pointer :: self

call mpas_geom_registry%init()
call mpas_geom_registry%add(c_key_self)
call mpas_geom_registry%get(c_key_self, self)

call geo_setup(self, c_conf)

end subroutine c_mpas_geo_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_clone(c_key_self, c_key_other) bind(c,name='mpas_geo_clone_f90')
use iso_c_binding
use mpas_geom_mod
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(mpas_geom), pointer :: self, other

call mpas_geom_registry%add(c_key_other)
call mpas_geom_registry%get(c_key_other, other)
call mpas_geom_registry%get(c_key_self, self)

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
