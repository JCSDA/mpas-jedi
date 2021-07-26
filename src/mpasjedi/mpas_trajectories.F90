! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_trajectories

use fckit_log_module, only: fckit_log

!MPAS-Model
use mpas_derived_types
use mpas_pool_routines

!mpas-jedi
use mpas_fields_mod, only: copy_pool, mpas_fields

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
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_traj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include <oops/util/linkedList_c.f>

! ------------------------------------------------------------------------------

subroutine set_traj(self,state)
implicit none
type(mpas_trajectory), intent(inout) :: self
type(mpas_fields),     intent(in)    :: state

 call fckit_log%info ('===> set_traj(self) in mpas_trajectories.F90')

 self%nf = state%nf
 allocate(self%fldnames( self%nf ))
 self%fldnames = state%fldnames
 call copy_pool(state % subFields, self % subFields)

 call fckit_log%info ('===> DONE set_traj(self) in mpas_trajectories.F90')

end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(mpas_trajectory), intent(inout) :: self

  call fckit_log%info ('===> delete_traj(self) in mpas_trajectories.F90')

end subroutine delete_traj

! ------------------------------------------------------------------------------

end module mpas_trajectories
