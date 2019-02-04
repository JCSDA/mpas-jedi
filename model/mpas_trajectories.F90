! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_trajectories

use mpas_derived_types

implicit none
private

public :: mpas_trajectory, set_traj, delete_traj
public :: mpas_traj_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold the model trajectory
type :: mpas_trajectory
  !integer :: nothing
  integer :: nf                                  ! Number of variables in fld
  character(len=22), allocatable  :: fldnames(:) ! Variable identifiers
  type (mpas_pool_type), pointer  :: subFields   !---> state variables (to be analyzed)
end type mpas_trajectory

#define LISTED_TYPE mpas_trajectory

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_traj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine set_traj(self,state)
use mpas_state_utils_mod
use mpas_field_utils_mod, only: copy_pool
use mpas_pool_routines
implicit none
type(mpas_trajectory), intent(inout) :: self
type(mpas_state)     , intent(in   ) :: state

 write(*,*) '===> set_traj(self) in mpas_trajectories.F90'

 self%nf = state%nf
 allocate(self%fldnames( self%nf ))
 self%fldnames = state%fldnames

 call copy_pool(state % subFields, self % subFields)

 write(*,*) '===> DONE set_traj(self) in mpas_trajectories.F90'

end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(mpas_trajectory), intent(inout) :: self

 write(*,*) '===> delete_traj(self) in mpas_trajectories.F90'

end subroutine delete_traj

! ------------------------------------------------------------------------------

end module mpas_trajectories
