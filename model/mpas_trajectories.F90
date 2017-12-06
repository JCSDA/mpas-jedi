! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_trajectories

implicit none
private

public :: mpas_trajectory, set_traj, delete_traj
public :: mpas_traj_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold the model trajectory
type :: mpas_trajectory
  integer :: nothing
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

subroutine set_traj(self)
implicit none
type(mpas_trajectory), intent(inout) :: self

end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(mpas_trajectory), intent(inout) :: self

end subroutine delete_traj

! ------------------------------------------------------------------------------

end module mpas_trajectories
