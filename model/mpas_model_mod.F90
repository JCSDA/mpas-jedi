! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_model_mod

use iso_c_binding
use config_mod
use duration_mod
use mpas_geom_mod
use mpas_fields_mod
use mpas_trajectories

implicit none
private
public :: mpas_model, & 
        & model_setup, model_delete, &
        & model_prepare_integration, model_prepare_integration_tl, model_prepare_integration_ad, &
        & model_propagate, model_propagate_tl, model_propagate_ad, &
        & model_prop_traj, model_wipe_traj, &
        & mpas_model_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: mpas_model
  integer :: nothing
end type mpas_model

#define LISTED_TYPE mpas_model

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine model_setup(model, geom, c_conf)
implicit none
type(c_ptr), intent(in) :: c_conf !< pointer to object of class Config
type(mpas_model)       :: model  ! should I put intent on these?
type(mpas_geom)        :: geom

end subroutine model_setup

! ------------------------------------------------------------------------------

subroutine model_delete(self)
implicit none
type(mpas_model) :: self

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_prepare_integration(self, flds)
implicit none
type(mpas_model) :: self
type(mpas_field) :: flds

end subroutine model_prepare_integration

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_ad(self, flds)
implicit none
type(mpas_model) :: self
type(mpas_field) :: flds

end subroutine model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_tl(self, flds)
implicit none
type(mpas_model) :: self
type(mpas_field) :: flds

end subroutine model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine model_propagate(self, flds)
implicit none
type(mpas_model) :: self
type(mpas_field) :: flds

end subroutine model_propagate

! ------------------------------------------------------------------------------

subroutine model_propagate_ad(self, flds, traj)
implicit none

type(mpas_model)      :: self
type(mpas_field)      :: flds
type(mpas_trajectory) :: traj

end subroutine model_propagate_ad

! ------------------------------------------------------------------------------

subroutine model_propagate_tl(self, flds, traj)
implicit none
type(mpas_model)      :: self
type(mpas_field)      :: flds
type(mpas_trajectory) :: traj

end subroutine model_propagate_tl

! ------------------------------------------------------------------------------

subroutine model_prop_traj(self, flds, traj)
implicit none
type(mpas_model)      :: self
type(mpas_field)      :: flds
type(mpas_trajectory) :: traj

end subroutine model_prop_traj

! ------------------------------------------------------------------------------

subroutine model_wipe_traj(traj)
implicit none
type(mpas_trajectory) :: traj

end subroutine model_wipe_traj

! ------------------------------------------------------------------------------

end module mpas_model_mod
