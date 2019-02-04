! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_mpas_varchange_setup(c_key_self, c_key_state_bg, c_key_state_fg, &
           c_key_geom, c_conf) bind (c,name='mpas_varchange_setup_f90')

use iso_c_binding
use mpas_varchange_mod
use mpas_state_utils_mod
use mpas_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_self     !< Change variable structure
integer(c_int), intent(in)    :: c_key_state_bg !< Background key
integer(c_int), intent(in)    :: c_key_state_fg !< First guess key
integer(c_int), intent(in)    :: c_key_geom     !< Geom key
type(c_ptr),    intent(in)    :: c_conf         !< Configuration

type(mpas_varchange), pointer :: self
type(mpas_state), pointer :: bg
type(mpas_state), pointer :: fg
type(mpas_geom), pointer :: geom

call mpas_varchange_registry%init()
call mpas_varchange_registry%add(c_key_self)
call mpas_varchange_registry%get(c_key_self, self)

call mpas_state_registry%get(c_key_state_bg,bg)
call mpas_state_registry%get(c_key_state_fg,fg)

call mpas_geom_registry%get(c_key_geom,geom)

call mpas_varchange_setup(self, bg, fg, geom, c_conf)

end subroutine c_mpas_varchange_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_varchange_delete(c_key_self) &
           bind (c,name='mpas_varchange_delete_f90')

use iso_c_binding
use mpas_varchange_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(mpas_varchange), pointer :: self

call mpas_varchange_registry%get(c_key_self,self)
call mpas_varchange_delete(self)
call mpas_varchange_registry%remove(c_key_self)

end subroutine c_mpas_varchange_delete

! ------------------------------------------------------------------------------

subroutine c_mpas_varchange_multiply(c_key_self, c_key_in, c_key_out) &
           bind (c,name='mpas_varchange_multiply_f90')

use iso_c_binding
use mpas_varchange_mod
use mpas_increment_utils_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_varchange), pointer :: self
type(mpas_increment), pointer :: xin
type(mpas_increment), pointer :: xout

call mpas_varchange_registry%get(c_key_self,self)
call mpas_increment_registry%get(c_key_in,xin)
call mpas_increment_registry%get(c_key_out,xout)

call mpas_varchange_multiply(self,xin,xout)

end subroutine c_mpas_varchange_multiply

! ----------------------------------------------------------------------------

subroutine c_mpas_varchange_multiplyadjoint(c_key_self, c_key_in, &
           c_key_out) bind (c,name='mpas_varchange_multiplyadjoint_f90')

use iso_c_binding
use mpas_varchange_mod
use mpas_increment_utils_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_varchange), pointer :: self
type(mpas_increment), pointer :: xin
type(mpas_increment), pointer :: xout

call mpas_varchange_registry%get(c_key_self,self)
call mpas_increment_registry%get(c_key_in,xin)
call mpas_increment_registry%get(c_key_out,xout)

call mpas_varchange_multiplyadjoint(self,xin,xout)

end subroutine c_mpas_varchange_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_mpas_varchange_multiplyinverse(c_key_self, c_key_in, c_key_out) &
           bind (c,name='mpas_varchange_multiplyinverse_f90')

use iso_c_binding
use mpas_varchange_mod
use mpas_increment_utils_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_varchange), pointer :: self
type(mpas_increment), pointer :: xin
type(mpas_increment), pointer :: xout

call mpas_varchange_registry%get(c_key_self,self)
call mpas_increment_registry%get(c_key_in,xin)
call mpas_increment_registry%get(c_key_out,xout)

call mpas_varchange_multiplyinverse(self,xin,xout)

end subroutine c_mpas_varchange_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_mpas_varchange_multiplyinverseadjoint(c_key_self, c_key_in, &
      c_key_out) bind (c,name='mpas_varchange_multiplyinverseadjoint_f90')

use iso_c_binding
use mpas_varchange_mod
use mpas_increment_utils_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpas_varchange), pointer :: self
type(mpas_increment), pointer :: xin
type(mpas_increment), pointer :: xout

call mpas_varchange_registry%get(c_key_self,self)
call mpas_increment_registry%get(c_key_in,xin)
call mpas_increment_registry%get(c_key_out,xout)

call mpas_varchange_multiplyinverseadjoint(self,xin,xout)

end subroutine c_mpas_varchange_multiplyinverseadjoint

! ----------------------------------------------------------------------------
