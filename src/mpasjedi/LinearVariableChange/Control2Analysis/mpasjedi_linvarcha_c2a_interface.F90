! (C) Copyright 2018-2021 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpasjedi_linvarcha_c2a_interface_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration

use oops_variables_mod

!use kinds
use mpas_fields_mod, only: mpas_fields, mpas_fields_registry
use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
use mpasjedi_linvarcha_c2a_mod !, only: mpasjedi_linvarcha_c2a

implicit none

private
public :: mpasjedi_linvarcha_c2a_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE mpasjedi_linvarcha_c2a

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpasjedi_linvarcha_c2a_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_create( &
           c_key_self, c_key_geom, c_key_bg, c_key_fg, c_conf, c_vars) &
           bind (c,name='mpasjedi_linvarcha_c2a_create_f90')

!implicit none
integer(c_int), intent(inout)  :: c_key_self     !< Change variable structure
integer(c_int), intent(in)     :: c_key_geom     !< Geom key
integer(c_int), intent(in)     :: c_key_bg       !< Background key
integer(c_int), intent(in)     :: c_key_fg       !< First guess key
type(c_ptr), value, intent(in) :: c_conf         !< Configuration
type(c_ptr), value, intent(in) :: c_vars         !< List of input variables

type(mpasjedi_linvarcha_c2a), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: bg
type(mpas_fields), pointer :: fg
type(fckit_configuration) :: f_conf
type(oops_variables) :: vars

call mpasjedi_linvarcha_c2a_registry%init()
call mpasjedi_linvarcha_c2a_registry%add(c_key_self)
call mpasjedi_linvarcha_c2a_registry%get(c_key_self, self)

call mpas_geom_registry%get(c_key_geom,geom)
call mpas_fields_registry%get(c_key_bg,bg)
call mpas_fields_registry%get(c_key_fg,fg)

f_conf = fckit_configuration(c_conf)
vars = oops_variables(c_vars)

call self%create(geom, bg, fg, f_conf, vars)

end subroutine c_mpasjedi_linvarcha_c2a_create

! ------------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_delete(c_key_self) &
           bind (c,name='mpasjedi_linvarcha_c2a_delete_f90')

!implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(mpasjedi_linvarcha_c2a), pointer :: self

call mpasjedi_linvarcha_c2a_registry%get(c_key_self,self)
call self%delete()
call mpasjedi_linvarcha_c2a_registry%remove(c_key_self)

end subroutine c_mpasjedi_linvarcha_c2a_delete

! ------------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_multiply( &
           c_key_self, c_key_geom, c_key_in, c_key_out) &
           bind (c,name='mpasjedi_linvarcha_c2a_multiply_f90')

!implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpasjedi_linvarcha_c2a), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: xin
type(mpas_fields), pointer :: xout

call mpasjedi_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_in,xin)
call mpas_fields_registry%get(c_key_out,xout)

call self%multiply(geom,xin,xout)

end subroutine c_mpasjedi_linvarcha_c2a_multiply

! ----------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_multiplyadjoint( &
           c_key_self, c_key_geom, c_key_in, &
           c_key_out) bind (c,name='mpasjedi_linvarcha_c2a_multiplyadjoint_f90')

!implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpasjedi_linvarcha_c2a), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: xin
type(mpas_fields), pointer :: xout

call mpasjedi_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_in,xin)
call mpas_fields_registry%get(c_key_out,xout)

call self%multiplyadjoint(geom,xin,xout)

end subroutine c_mpasjedi_linvarcha_c2a_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_multiplyinverse( &
           c_key_self, c_key_geom, c_key_in, c_key_out) &
           bind (c,name='mpasjedi_linvarcha_c2a_multiplyinverse_f90')

!implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpasjedi_linvarcha_c2a), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: xin
type(mpas_fields), pointer :: xout

call mpasjedi_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_in,xin)
call mpas_fields_registry%get(c_key_out,xout)

call self%multiplyinverse(geom,xin,xout)

end subroutine c_mpasjedi_linvarcha_c2a_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_mpasjedi_linvarcha_c2a_multiplyinverseadjoint( &
           c_key_self, c_key_geom, c_key_in, &
           c_key_out) bind (c,name='mpasjedi_linvarcha_c2a_multiplyinverseadjoint_f90')

!implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(mpasjedi_linvarcha_c2a), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: xin
type(mpas_fields), pointer :: xout

call mpasjedi_linvarcha_c2a_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_in,xin)
call mpas_fields_registry%get(c_key_out,xout)

call self%multiplyinverseadjoint(geom,xin,xout)

end subroutine c_mpasjedi_linvarcha_c2a_multiplyinverseadjoint

! ----------------------------------------------------------------------------

end module mpasjedi_linvarcha_c2a_interface_mod
