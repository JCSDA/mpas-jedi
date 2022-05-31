! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module mpasjedi_vc_model2geovars_interface_mod

use iso_c_binding

use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use mpas_geom_mod, only: mpas_geom, mpas_geom_registry
use mpas_fields_mod, only: mpas_fields, mpas_fields_registry
use mpasjedi_vc_model2geovars_mod, only: mpasjedi_vc_model2geovars

implicit none

private
public :: mpasjedi_vc_model2geovars_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE mpasjedi_vc_model2geovars

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpasjedi_vc_model2geovars_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_vc_model2geovars_create(c_key_self, c_key_geom, c_conf) &
           bind (c, name='mpasjedi_vc_model2geovars_create_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_key_geom
type(c_ptr), value, intent(in) :: c_conf

type(mpasjedi_vc_model2geovars), pointer :: self
type(mpas_geom), pointer :: geom
type(fckit_configuration) :: conf

! Linked list
! -----------
call mpasjedi_vc_model2geovars_registry%init()
call mpasjedi_vc_model2geovars_registry%add(c_key_self)
call mpasjedi_vc_model2geovars_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom,geom)

! APIs
! ----
conf = fckit_configuration(c_conf)

! Implementation
! --------------
call self%create(geom, conf)

end subroutine c_mpasjedi_vc_model2geovars_create

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_vc_model2geovars_delete(c_key_self) &
           bind (c, name='mpasjedi_vc_model2geovars_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(mpasjedi_vc_model2geovars), pointer :: self

! Linked list
! -----------
call mpasjedi_vc_model2geovars_registry%get(c_key_self,self)

! Implementation
! --------------
!call self%delete()

! Linked list
! -----------
call mpasjedi_vc_model2geovars_registry%remove(c_key_self)

end subroutine c_mpasjedi_vc_model2geovars_delete

! --------------------------------------------------------------------------------------------------

subroutine c_mpasjedi_vc_model2geovars_changevar(c_key_self, c_key_geom, c_key_xm, c_key_xg) &
           bind (c, name='mpasjedi_vc_model2geovars_changevar_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xm
integer(c_int), intent(in) :: c_key_xg

type(mpasjedi_vc_model2geovars), pointer :: self
type(mpas_geom), pointer :: geom
type(mpas_fields), pointer :: xm
type(mpas_fields), pointer :: xg

! Linked list
! -----------
call mpasjedi_vc_model2geovars_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_xm,xm)
call mpas_fields_registry%get(c_key_xg,xg)
call mpas_geom_registry%get(c_key_geom,geom)

! Implementation
! --------------
call self%changevar(geom, xm, xg)

end subroutine c_mpasjedi_vc_model2geovars_changevar

! --------------------------------------------------------------------------------------------------

end module mpasjedi_vc_model2geovars_interface_mod
