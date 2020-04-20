! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module mpasjedi_getvalues_interface_mod

! Intrinsic
use iso_c_binding
use kinds,               only: kind_real

! oops dependencies
use datetime_mod
use duration_mod
use oops_variables_mod

! ufo dependencies
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

! self dependency
use mpasjedi_getvalues_mod, only: mpasjedi_getvalues, mpas_getvalues_registry

! mpas dependencies
use mpas_geom_mod,            only: mpas_geom_registry, mpas_geom
use mpas_field_utils_mod
! use mpas_state_interface_mod, only: mpas_state_registry

implicit none
private

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine mpas_getvalues_create_c(c_key_self, c_key_geom, c_key_locs) &
           bind (c, name='mpas_getvalues_create_f90')
implicit none
integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
integer(c_int),     intent(in)    :: c_key_locs      !< Key to observation locations

type(mpasjedi_getvalues),  pointer :: self
type(mpas_geom),           pointer :: geom
type(ufo_locs),            pointer :: locs

! Create object
call mpas_getvalues_registry%init()
call mpas_getvalues_registry%add(c_key_self)
call mpas_getvalues_registry%get(c_key_self, self)

! Others
call mpas_geom_registry%get(c_key_geom, geom)
call ufo_locs_registry%get(c_key_locs, locs)

! Call method
call self%create(geom, locs)

end subroutine mpas_getvalues_create_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_getvalues_delete_c(c_key_self) bind (c, name='mpas_getvalues_delete_f90')

integer(c_int), intent(inout) :: c_key_self !< Key to self

type(mpasjedi_getvalues), pointer :: self

! Get object
call mpas_getvalues_registry%get(c_key_self, self)

! Call method
call self%delete()

! Remove object
call mpas_getvalues_registry%remove(c_key_self)

end subroutine mpas_getvalues_delete_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_getvalues_fill_geovals_c(c_key_self, c_key_geom, c_key_state, c_t1, c_t2, &
                                            c_key_locs, c_key_geovals) &
           bind (c, name='mpas_getvalues_fill_geovals_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
integer(c_int),     intent(in) :: c_key_locs
integer(c_int),     intent(in) :: c_key_geovals

type(mpasjedi_getvalues), pointer :: self
type(mpas_geom),          pointer :: geom
type(mpas_field),         pointer :: fields
type(datetime)                    :: t1
type(datetime)                    :: t2
type(ufo_locs),           pointer :: locs
type(ufo_geovals),        pointer :: geovals

! Get objects
call mpas_getvalues_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_field_registry%get(c_key_state, fields)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
call ufo_locs_registry%get(c_key_locs, locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, fields, t1, t2, locs, geovals)

end subroutine mpas_getvalues_fill_geovals_c

! --------------------------------------------------------------------------------------------------

end module mpasjedi_getvalues_interface_mod

