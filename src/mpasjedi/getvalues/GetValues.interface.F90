! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module mpasjedi_getvalues_interface_mod

! Intrinsic
use iso_c_binding
use kinds,               only: kind_real
use fckit_configuration_module, only: fckit_configuration

! oops dependencies
use datetime_mod

! ufo dependencies
use ufo_locations_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

! self dependency
use mpasjedi_getvalues_mod, only: mpasjedi_getvalues, mpas_getvalues_registry

! mpas dependencies
use mpas_geom_mod,            only: mpas_geom_registry, mpas_geom
use mpas_fields_mod
! use mpas_state_interface_mod, only: mpas_state_registry

implicit none
private

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine mpas_getvalues_create_c(c_key_self, c_key_geom, c_locs, c_conf) &
           bind (c, name='mpas_getvalues_create_f90')
implicit none
integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
type(c_ptr), value, intent(in)    :: c_locs          !< Observation locations
type(c_ptr), value, intent(in)    :: c_conf          !< Key to configuration

type(mpasjedi_getvalues),  pointer :: self
type(mpas_geom),           pointer :: geom
type(ufo_locations)                :: locs
type(fckit_configuration)          :: f_conf

! Create object
call mpas_getvalues_registry%init()
call mpas_getvalues_registry%add(c_key_self)
call mpas_getvalues_registry%get(c_key_self, self)

! Others
call mpas_geom_registry%get(c_key_geom, geom)
locs = ufo_locations(c_locs)
f_conf = fckit_configuration(c_conf)

! Call method
call self%create(geom, locs, f_conf)

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
                                         c_locs, c_key_geovals) &
           bind (c, name='mpas_getvalues_fill_geovals_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(mpasjedi_getvalues), pointer :: self
type(mpas_geom),          pointer :: geom
type(mpas_fields),        pointer :: fields
type(datetime)                    :: t1
type(datetime)                    :: t2
type(ufo_locations)               :: locs
type(ufo_geovals),        pointer :: geovals

! Get objects
call mpas_getvalues_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_state, fields)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, fields, t1, t2, locs, geovals)

end subroutine mpas_getvalues_fill_geovals_c

! --------------------------------------------------------------------------------------------------

end module mpasjedi_getvalues_interface_mod

