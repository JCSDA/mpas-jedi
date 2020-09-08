! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module mpas_linvarcha_c2a_interface_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding

use kinds
use mpas_field_utils_mod
use mpas_geom_mod
use mpas_linvarcha_c2a_mod

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_setup(c_key_self, c_key_state_bg, c_key_state_fg, &
           c_key_geom, c_conf) bind (c,name='mpas_linvarcha_c2a_setup_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self     !< Change variable structure
integer(c_int), intent(in)     :: c_key_state_bg !< Background key
integer(c_int), intent(in)     :: c_key_state_fg !< First guess key
integer(c_int), intent(in)     :: c_key_geom     !< Geom key
type(c_ptr), value, intent(in) :: c_conf         !< Configuration

type(mpas_linvarcha_c2a), pointer :: self
type(mpas_field), pointer :: bg
type(mpas_field), pointer :: fg
type(mpas_geom), pointer :: geom
type(fckit_configuration) :: f_conf

call mpas_linvarcha_c2a_registry%init()
call mpas_linvarcha_c2a_registry%add(c_key_self)
call mpas_linvarcha_c2a_registry%get(c_key_self, self)

call mpas_field_registry%get(c_key_state_bg,bg)
call mpas_field_registry%get(c_key_state_fg,fg)

call mpas_geom_registry%get(c_key_geom,geom)

f_conf = fckit_configuration(c_conf)
call mpas_linvarcha_c2a_setup(self, bg, fg, geom, f_conf)

end subroutine c_mpas_linvarcha_c2a_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_delete(c_key_self) &
           bind (c,name='mpas_linvarcha_c2a_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(mpas_linvarcha_c2a), pointer :: self

call mpas_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_linvarcha_c2a_delete(self)
call mpas_linvarcha_c2a_registry%remove(c_key_self)

end subroutine c_mpas_linvarcha_c2a_delete

! ------------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_multiply(c_key_self, c_key_in, c_key_out) &
           bind (c,name='mpas_linvarcha_c2a_multiply_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_linvarcha_c2a), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)

call mpas_linvarcha_c2a_multiply(self,xin,xout)

end subroutine c_mpas_linvarcha_c2a_multiply

! ----------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_multiplyadjoint(c_key_self, c_key_in, &
           c_key_out) bind (c,name='mpas_linvarcha_c2a_multiplyadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_linvarcha_c2a), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)

call mpas_linvarcha_c2a_multiplyadjoint(self,xin,xout)

end subroutine c_mpas_linvarcha_c2a_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_multiplyinverse(c_key_self, c_key_in, c_key_out) &
           bind (c,name='mpas_linvarcha_c2a_multiplyinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_linvarcha_c2a), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)

call mpas_linvarcha_c2a_multiplyinverse(self,xin,xout)

end subroutine c_mpas_linvarcha_c2a_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_mpas_linvarcha_c2a_multiplyinverseadjoint(c_key_self, c_key_in, &
      c_key_out) bind (c,name='mpas_linvarcha_c2a_multiplyinverseadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_linvarcha_c2a), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)

call mpas_linvarcha_c2a_multiplyinverseadjoint(self,xin,xout)

end subroutine c_mpas_linvarcha_c2a_multiplyinverseadjoint

! ----------------------------------------------------------------------------

end module mpas_linvarcha_c2a_interface_mod
