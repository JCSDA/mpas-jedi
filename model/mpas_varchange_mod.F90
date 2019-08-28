! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_varchange_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds

use mpas_state_utils_mod, only: mpas_state
use mpas_increment_utils_mod, only: mpas_increment
use mpas_geom_mod,   only: mpas_geom

implicit none

!> Fortran derived type to hold configuration data for the B mat variable change
type :: mpas_varchange
 integer :: tablesize
 real(8), allocatable :: estblx(:)
end type mpas_varchange

#define LISTED_TYPE mpas_varchange

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_varchange_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine mpas_varchange_setup(self, bg, fg, geom, f_conf)

implicit none
type(mpas_varchange),      intent(inout) :: self    !< Change variable structure
type(mpas_state), target,  intent(in)    :: bg
type(mpas_state), target,  intent(in)    :: fg
type(mpas_geom),           intent(in)    :: geom
type(fckit_configuration), intent(in)    :: f_conf  !< Configuration

end subroutine mpas_varchange_setup

! ------------------------------------------------------------------------------

subroutine mpas_varchange_delete(self)

implicit none
type(mpas_varchange), intent(inout) :: self

end subroutine mpas_varchange_delete

! ------------------------------------------------------------------------------

subroutine mpas_varchange_multiply(self,xctl,xmod)

implicit none
type(mpas_varchange), intent(inout) :: self
type(mpas_increment), intent(inout) :: xctl
type(mpas_increment), intent(inout) :: xmod

end subroutine mpas_varchange_multiply

! ------------------------------------------------------------------------------

subroutine mpas_varchange_multiplyadjoint(self,xmod,xctl)

implicit none
type(mpas_varchange), intent(inout) :: self
type(mpas_increment), intent(inout) :: xmod
type(mpas_increment), intent(inout) :: xctl

!Adjoint of analysis (control) to model variables

end subroutine mpas_varchange_multiplyadjoint

! ------------------------------------------------------------------------------

subroutine mpas_varchange_multiplyinverse(self,xinc,xctr)

implicit none
type(mpas_varchange), intent(inout) :: self
type(mpas_increment), intent(inout) :: xinc
type(mpas_increment), intent(inout) :: xctr

!> Not implemented

end subroutine mpas_varchange_multiplyinverse

! ------------------------------------------------------------------------------

subroutine mpas_varchange_multiplyinverseadjoint(self,xinc,xctr)

implicit none
type(mpas_varchange), intent(inout) :: self
type(mpas_increment), intent(inout) :: xinc
type(mpas_increment), intent(inout) :: xctr

!> Not implemented

end subroutine mpas_varchange_multiplyinverseadjoint

! ------------------------------------------------------------------------------

end module mpas_varchange_mod
