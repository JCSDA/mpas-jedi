! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_lvc_model2geovars_interface_mod

use iso_c_binding

use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
use mpas_fields_mod, only: mpas_fields, mpas_fields_registry
use mpasjedi_lvc_model2geovars_mod, only: mpasjedi_lvc_model2geovars

implicit none

private
public :: mpasjedi_lvc_model2geovars_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE mpasjedi_lvc_model2geovars

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpasjedi_lvc_model2geovars_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_lvc_model2geovars_create( &
  c_key_self, c_key_geom, c_key_bg, c_key_fg, c_conf) &
  bind (c,name='mpasjedi_lvc_model2geovars_create_f90')

integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in)    :: c_key_geom
integer(c_int), intent(in)    :: c_key_bg
integer(c_int), intent(in)    :: c_key_fg
type(c_ptr), value, intent(in) :: c_conf

type(mpasjedi_lvc_model2geovars), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: bg
type(mpas_fields), pointer :: fg
type(fckit_configuration) :: conf

! Linked list
! -----------
call mpasjedi_lvc_model2geovars_registry%init()
call mpasjedi_lvc_model2geovars_registry%add(c_key_self)
call mpasjedi_lvc_model2geovars_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_bg, bg)
call mpas_fields_registry%get(c_key_fg, fg)

! APIs
! ----
conf = fckit_configuration(c_conf)

! Implementation
! --------------
call self%create(geom, bg, fg, conf)

end subroutine c_mpasjedi_lvc_model2geovars_create

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_lvc_model2geovars_delete(c_key_self) &
  bind (c,name='mpasjedi_lvc_model2geovars_delete_f90')

integer(c_int), intent(inout) :: c_key_self

type(mpasjedi_lvc_model2geovars), pointer :: self

! Linked list
! -----------
call mpasjedi_lvc_model2geovars_registry%get(c_key_self, self)

! Implementation
! --------------
call self%delete()

! Linked list
! -----------
call mpasjedi_lvc_model2geovars_registry%remove(c_key_self)

end subroutine c_mpasjedi_lvc_model2geovars_delete

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_lvc_model2geovars_multiply( &
  c_key_self, c_key_geom, c_key_dxm, c_key_dxg) &
  bind (c,name='mpasjedi_lvc_model2geovars_multiply_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_dxm
integer(c_int), intent(in) :: c_key_dxg

type(mpasjedi_lvc_model2geovars), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: dxm
type(mpas_fields), pointer :: dxg

! Linked list
! -----------
call mpasjedi_lvc_model2geovars_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom,geom)
call mpas_fields_registry%get(c_key_dxm,dxm)
call mpas_fields_registry%get(c_key_dxg,dxg)

! Implementation
! --------------
call self%multiply(geom, dxm, dxg)

end subroutine c_mpasjedi_lvc_model2geovars_multiply

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_lvc_model2geovars_multiplyadjoint( &
  c_key_self, c_key_geom, c_key_dxg, c_key_dxm) &
  bind (c,name='mpasjedi_lvc_model2geovars_multiplyadjoint_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_dxg
integer(c_int), intent(in) :: c_key_dxm

type(mpasjedi_lvc_model2geovars), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: dxg
type(mpas_fields), pointer :: dxm

! Linked list
! -----------
call mpasjedi_lvc_model2geovars_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%get(c_key_dxg, dxg)
call mpas_fields_registry%get(c_key_dxm, dxm)

! Implementation
! --------------
call self%multiplyadjoint(geom, dxg, dxm)

end subroutine c_mpasjedi_lvc_model2geovars_multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module mpasjedi_lvc_model2geovars_interface_mod
