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

use mpas_derived_types
use mpas_framework
use mpas_kind_types
use mpas_subdriver
use atm_core
use mpi 
use mpas_stream_manager
use mpas2da_mod

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
   type (domain_type), pointer :: domain 
   type (core_type), pointer :: corelist
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

subroutine model_setup(self, geom, c_conf)
implicit none
type(c_ptr), intent(in) :: c_conf !< pointer to object of class Config
type(mpas_model), intent(inout) :: self  ! should I put intent on these?
type(mpas_geom)        :: geom

character(len=20) :: ststep
type(duration) :: dtstep

write(0,*)'geom % nCells: ',geom % nCells
!self % corelist => geom % corelist
!self % domain => geom % domain

write(*,*)'===> model_setup'

ststep = config_get_string(c_conf,len(ststep),"tstep")
write(0,*)'ststep, dtstep: ', ststep

end subroutine model_setup

! ------------------------------------------------------------------------------

subroutine model_delete(self)
implicit none
type(mpas_model) :: self

write(*,*)'===> model_delete'

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_prepare_integration(self, flds)
implicit none
type(mpas_model) :: self
type(mpas_field) :: flds

write(*,*)'===> model_prepare_integration'

!call da_copy_sub2all_fields(self % domain, self % flds)

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
type (mpas_pool_type), pointer :: state

real (kind=RKIND) :: dt
integer :: itimestep, ii

write(*,*)'===> model_propagate'

!dt = 1200.
!itimestep = 1
!call mpas_pool_get_subpool(self % domain % blocklist % structs,'state',state)
!
!do ii=1, 1
!   
!   call atm_do_timestep(self % domain, dt, ii)
!   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
!   call mpas_pool_shift_time_levels(state)
!   itimestep = itimestep + 1
!   call mpas_advance_clock(clock)
!
!end do
!
!call da_copy_all2sub_fields(self % domain, flds % subFields)

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
