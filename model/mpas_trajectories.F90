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
  !type (mpas_pool_type)  :: subFields   !---> state variables (to be analyzed)
  !type (mpas_pool_type)  :: auxFields   !---> auxiliary variables, such as pressure, t2m, u10, v10, Tsfc
  type (mpas_pool_type), pointer  :: subFields   !---> state variables (to be analyzed)
  type (mpas_pool_type), pointer  :: auxFields   !---> auxiliary variables, such as pressure, t2m, u10, v10, Tsfc
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

subroutine set_traj(self,flds)
use mpas_fields_mod
use mpas_pool_routines
implicit none
type(mpas_trajectory), intent(inout) :: self
type(mpas_field)     , intent(in   ) :: flds

 write(*,*) '===> set_traj(self) in mpas_trajectories.F90'

 self%nf = flds%nf
 allocate(self%fldnames( self%nf ))
 self%fldnames = flds%fldnames

 write(*,*) ' associated(self%subFields) = ', associated(self%subFields) 
 if( associated(self%subFields) ) then
   call mpas_pool_empty_pool(self % subFields)
   call mpas_pool_destroy_pool(self % subFields)
   call mpas_pool_empty_pool(self % auxFields)
   call mpas_pool_destroy_pool(self % auxFields)
 endif

 call mpas_pool_create_pool(self % subFields,self % nf)
 call mpas_pool_clone_pool(flds % subFields, self % subFields)
 call mpas_pool_create_pool(self % auxFields) !???,self % nf)
 call mpas_pool_clone_pool(flds % subFields, self % auxFields)

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
