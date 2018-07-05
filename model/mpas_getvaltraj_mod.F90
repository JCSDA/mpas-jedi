
!> Fortran module handling interpolation trajectory for the MPAS model

module mpas_getvaltraj_mod

!General JEDI uses
use kinds
use iso_c_binding
use type_bump, only: bump_type
use mpas_pool_routines

implicit none
private

public mpas_getvaltraj
public mpas_getvaltraj_registry
public c_mpas_getvaltraj_setup, c_mpas_getvaltraj_delete

type :: mpas_getvaltraj
 integer :: nobs, ngrid
 type (mpas_pool_type), pointer :: pool_traj
 type(bump_type) :: bump
 integer :: nsize = 3 !< size of pool, currently for theta, index_qc, pressure
 logical :: lalloc = .false.
end type mpas_getvaltraj

#define LISTED_TYPE mpas_getvaltraj

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_getvaltraj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_mpas_getvaltraj_setup(c_key_self) bind(c,name='mpas_getvaltraj_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_getvaltraj), pointer :: self
integer :: nsize  !< size of trajectory pool

! Init, add and get key
! ---------------------
call mpas_getvaltraj_registry%init()
call mpas_getvaltraj_registry%add(c_key_self)
call mpas_getvaltraj_registry%get(c_key_self,self)

print*, 'dh: getvaltraj_setup', c_key_self

self%lalloc = .false.
self%nobs = 0
self%ngrid = 0

call mpas_pool_create_pool(self % pool_traj, self % nsize)

end subroutine c_mpas_getvaltraj_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_getvaltraj_delete(c_key_self) bind(c,name='mpas_getvaltraj_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_getvaltraj), pointer :: self

! Get key
call mpas_getvaltraj_registry%get(c_key_self, self)

if (self%lalloc) then
  self%nobs = 0
  self%ngrid = 0
  call mpas_pool_empty_pool( self % pool_traj )
  call mpas_pool_destroy_pool( self % pool_traj )
  call self%bump%dealloc
  self%lalloc = .false.
endif

! Remove key
call mpas_getvaltraj_registry%remove(c_key_self)

end subroutine c_mpas_getvaltraj_delete

! ------------------------------------------------------------------------------

end module mpas_getvaltraj_mod
