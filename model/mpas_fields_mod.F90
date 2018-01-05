! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_fields_mod

!use iso_c_binding , only: c_char, c_int
use iso_c_binding 
use config_mod
use datetime_mod
use mpas_geom_mod
use mpas_vars_mod
use mpas_kinds
use ufo_locs_mod
use ufo_geovals_mod

use mpas_derived_types
use mpas_framework
use mpas_kind_types
!use init_atm_core_interface
use mpas_subdriver
use atm_core
use mpas2da_mod

implicit none
private

public :: mpas_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp_tl, interp_ad
public :: mpas_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold MPAS fields
type :: mpas_field
  type (domain_type), pointer :: domain 
  type (core_type), pointer :: corelist
  !type (dm_info), pointer :: dminfo
  type (mpas_geom), pointer :: geom                                 !< Number of unstructured grid cells
  integer :: nf                           !< Number of variables in fld
  !integer :: sum_scalar                   !< Number of variables in fld
  !integer :: sum_aero                     !< Number of variables in fld
  !integer :: ns                           !< Number of surface fields (x1d [nCells])
  character(len=20), allocatable :: fldnames(:)      !< Variable identifiers
  type (mpas_pool_type), pointer :: subFields
end type mpas_field

#define LISTED_TYPE mpas_field

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_field_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)


    implicit none

    type(mpas_field), intent(inout) :: self
    type(mpas_geom),  intent(inout)    :: geom
    type(mpas_vars),  intent(in)    :: vars

    type (dm_info), pointer :: dminfo
    integer :: nsize

    ! from the namelist
    self % nf = vars % nv
    allocate(self % fldnames(self % nf))
    self % fldnames(:) = vars % fldnames(:)

    ! coming from mpas_subdriver
    call mpas_init(self % corelist, self % domain)

    ! update geom
    self % geom = geom
    ! This or create a subpool dminfo and clone
    geom % dminfo => self % domain % dminfo

    ! Create a subpool from allfields
    nsize = da_common_vars(self % domain % blocklist % allFields, self % fldnames)
    if ( self % nf .ne. nsize  ) then
       call abor1_ftn("mpas_fields:create: dimension mismatch ",self % nf, nsize)
    end if
    write(0,*)'-- Create a sub Pool from list of variable ',nsize
    call mpas_pool_create_pool(self % subFields,self % nf)
    call da_make_subpool(self % domain % blocklist % allFields, self % subFields, self % fldnames)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(mpas_field), intent(inout) :: self

   write(0,*)'--> deallocate subFields Pool'
   call mpas_pool_empty_pool(self % subFields)
   call mpas_pool_destroy_pool(self % subFields)
   write(0,*)'--> deallocate domain and core'
   call mpas_finalize(self % corelist, self % domain)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(mpas_field), intent(inout) :: self

  call da_zeros(self % subFields)

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine random(self)
implicit none
type(mpas_field), intent(inout) :: self

  !call da_random(self % subFields)

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

      !call mpas_pool_clone_pool(da_state, da_state_incr)

end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
character(len=StrKIND) :: kind_op

kind_op = 'add'
call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
character(len=StrKIND) :: kind_op

kind_op = 'schur'
call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
character(len=StrKIND) :: kind_op

kind_op = 'sub'
call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(mpas_field), intent(inout)  :: self
real(kind=kind_real), intent(in) :: zz

call da_self_mult(self % subFields, zz)

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(mpas_field), intent(inout)  :: self
real(kind=kind_real), intent(in) :: zz
type(mpas_field), intent(in)     :: rhs

call da_axpy(self % subFields, rhs % subFields, zz)

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)
implicit none
type(mpas_field), intent(in) :: fld1, fld2
real(kind=kind_real), intent(inout) :: zprod

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
character(len=StrKIND) :: kind_op

! GD: I don't see any difference than for self_add other than subFields can contain
! different variables than mpas_field and the resolution of incr can be different. 

if (self%geom%nCells==rhs%geom%nCells .and. self%geom%nVertLevels==rhs%geom%nVertLevels) then
   kind_op = 'add'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields)
else
   call abor1_ftn("mpas_fields:add_incr: dimension mismatch")
endif

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(mpas_field), intent(inout) :: lhs
type(mpas_field), intent(in)    :: x1
type(mpas_field), intent(in)    :: x2
character(len=StrKIND) :: kind_op

call zeros(lhs)
if (x1%geom%nCells==x2%geom%nCells .and. x1%geom%nVertLevels==x2%geom%nVertLevels) then
  if (lhs%geom%nCells==x1%geom%nCells .and. lhs%geom%nVertLevels==x1%geom%nVertLevels) then
     kind_op = 'sub'
     call da_operator(trim(kind_op), lhs % subFields, x1 % subFields, x2 % subFields)
  else
    call abor1_ftn("mpas_fields:diff_incr: dimension mismatch between the two variables.")
  endif
else
  call abor1_ftn("mpas_fields:diff_incr: states not at same resolution")
endif

return

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(mpas_field), intent(inout) :: fld
type(mpas_field), intent(in)    :: rhs

! FIXME: We just copy rhs to fld for now. Need an actual interpolation routine later. (SH)
if (fld%geom%nCells == rhs%geom%nCells .and.  fld%geom%nVertLevels == rhs%geom%nVertLevels) then
  call copy(fld, rhs)
else
  write(0,*) fld%geom%nCells, rhs%geom%nCells, fld%geom%nVertLevels, rhs%geom%nVertLevels
  call abor1_ftn("mpas_fields:field_resol: dimension mismatch")
endif

end subroutine change_resol

! -----------------------------------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)

!use mpas_derived_types
!use mpas_pool_routines
!use mpas_dmpar
!use mpas_abort, only : mpas_dmpar_global_abort
!use mpas_stream_manager

!!type (MPAS_Clock_type), pointer :: clock
implicit none
type(mpas_field), intent(inout) :: fld      !< Fields
type(c_ptr), intent(in)          :: c_conf   !< Configuration
type(datetime), intent(inout)    :: vdate    !< DateTime

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)
implicit none
type(mpas_field), intent(in) :: fld    !< Fields
type(c_ptr), intent(in)       :: c_conf !< Configuration
type(datetime), intent(in)    :: vdate  !< DateTime

!call mpas_dmpar_get_time(output_start_time)
!call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)
implicit none
type(mpas_field), intent(in) :: fld
integer, intent(in) :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
implicit none
type(mpas_field), intent(in)      :: fld
real(kind=kind_real), intent(out) :: prms

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, vars, gom)
implicit none
type(mpas_field), intent(in)    :: fld
type(ufo_locs), intent(in)  :: locs
type(mpas_vars), intent(in)     :: vars
type(ufo_geovals), intent(inout) :: gom

end subroutine interp_tl

! ------------------------------------------------------------------------------

subroutine interp_ad(fld, locs, vars, gom)
implicit none
type(mpas_field), intent(inout) :: fld
type(ufo_locs), intent(in)  :: locs
type(mpas_vars), intent(in)     :: vars
type(ufo_geovals), intent(inout) :: gom

end subroutine interp_ad

! ------------------------------------------------------------------------------

end module mpas_fields_mod
