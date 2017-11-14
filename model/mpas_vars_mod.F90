
!> Fortran module to handle variables for the MPAS-A model
module mpas_vars_mod

use fckit_log_module, only : log
use iso_c_binding
use config_mod

implicit none
private
public :: mpas_vars, mpas_vars_setup, mpas_vars_clone
public :: mpas_vars_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to represent MPAS-A model variables
type :: mpas_vars
  integer :: nv
  character(len=22), allocatable :: fldnames(:) !< Variable identifiers
  logical :: lbc
end type mpas_vars

#define LISTED_TYPE mpas_vars

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_vars_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine mpas_vars_setup(self, cvars)
implicit none
type(mpas_vars), intent(inout) :: self
character(len=22), intent(in) :: cvars(:)
integer :: jj

self%nv = size(cvars)
self%lbc = .false.

do jj=1,self%nv
  if (trim(cvars(jj))/="theta" .and. trim(cvars(jj))/="rho" .and. trim(cvars(jj))/="qv" .and. &
      trim(cvars(jj))/="u" .and. trim(cvars(jj))/="uReconstructZonal" .and. trim(cvars(jj))/="uReconstructMeridional") &
      call abor1_ftn ("mpas_vars_setup: unknown field")
enddo
allocate(self%fldnames(self%nv))
self%fldnames(:)=cvars(:)

end subroutine mpas_vars_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_vars_create(c_key_self, c_conf) bind(c,name='mpas_var_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(mpas_vars), pointer :: self
character(len=20) :: svar

call mpas_vars_registry%init()
call mpas_vars_registry%add(c_key_self)
call mpas_vars_registry%get(c_key_self, self)

svar = config_get_string(c_conf,len(svar),"variables")
call log%info('SVAR: '//svar)
select case (trim(svar))
case ("reconstructed_winds")
  self%nv = 5
  self%lbc = .false.
  allocate(self%fldnames(self%nv))
!  self%fldnames(:) = (/"theta","rho","qv","uReconstructZonal","uReconstructMeridional"/)
  self%fldnames(1) = "theta"
  self%fldnames(2) = "rho"
  self%fldnames(3) = "qv"
  self%fldnames(4) = "uReconstructZonal"
  self%fldnames(5) = "uReconstructMeridional"
case ("normal_speed")
  self%nv = 4
  self%lbc = .false.
  allocate(self%fldnames(self%nv))
  !self%fldnames(:) = (/"theta","rho","qv","u"/)
  self%fldnames(1) = "theta"
  self%fldnames(2) = "rho"
  self%fldnames(3) = "qv"
  self%fldnames(4) = "u"
case ("onevar")
  self%nv = 1
  self%lbc = .false.
  allocate(self%fldnames(self%nv))
  self%fldnames(1) = "theta"
case default
  call abor1_ftn("c_mpas_vars_create: undefined variables")
end select

return
end subroutine c_mpas_vars_create

! ------------------------------------------------------------------------------

subroutine c_mpas_vars_clone(c_key_self, c_key_other) bind(c,name='mpas_var_clone_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(mpas_vars), pointer :: self, other

call mpas_vars_registry%get(c_key_self, self)
call mpas_vars_registry%add(c_key_other)
call mpas_vars_registry%get(c_key_other, other)

call mpas_vars_clone(self, other)

end subroutine c_mpas_vars_clone

! ------------------------------------------------------------------------------

subroutine mpas_vars_clone(self, other)
implicit none
type(mpas_vars), intent(in)    :: self
type(mpas_vars), intent(inout) :: other

other%nv = self%nv
other%lbc = self%lbc

allocate(other%fldnames(other%nv))
other%fldnames(:)=self%fldnames(:)

end subroutine mpas_vars_clone

! ------------------------------------------------------------------------------

subroutine c_mpas_vars_delete(c_key_self) bind(c,name='mpas_var_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(mpas_vars), pointer :: self
call mpas_vars_registry%get(c_key_self, self)
deallocate(self%fldnames)
call mpas_vars_registry%remove(c_key_self)

return
end subroutine c_mpas_vars_delete

! ------------------------------------------------------------------------------

subroutine c_mpas_vars_info(c_key_self, c_nv, c_nl) bind(c,name='mpas_var_info_f90')
implicit none
integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_nv
integer(c_int), intent(inout) :: c_nl
type(mpas_vars), pointer :: self

call mpas_vars_registry%get(c_key_self, self)

c_nv = self%nv
c_nl = 0
if (self%lbc) c_nl = 1

!self%fldnames(:)

return
end subroutine c_mpas_vars_info

! ------------------------------------------------------------------------------

end module mpas_vars_mod
