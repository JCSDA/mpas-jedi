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
!use mpas_subdriver
!use mpas_core
!use mpas_subdriver


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
  integer :: sum_scalar                   !< Number of variables in fld
  integer :: sum_aero                     !< Number of variables in fld
  integer :: ns                           !< Number of surface fields (x1d [nCells])
  character(len=20), allocatable :: fldnames(:)      !< Variable identifiers
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
    type(mpas_geom),  intent(in)    :: geom
    type(mpas_vars),  intent(in)    :: vars

    type (core_type), pointer :: corelist => null()
    type (dm_info), pointer :: dminfo
    type (domain_type), pointer :: domain_ptr

 
    !--- it compiles --------------------------
    !allocate(corelist)
    !nullify(corelist % next)

    !allocate(corelist % domainlist)
    !nullify(corelist % domainlist % next)

    !domain_ptr => corelist % domainlist
    !domain_ptr % core => corelist

    !call mpas_allocate_domain(domain_ptr)
    ! --- end of it compiles ------------------

    !call mpas_init()

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(mpas_field), intent(inout) :: self

!call  mpas_finalize()

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(mpas_field), intent(inout) :: self

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine random(self)
implicit none
type(mpas_field), intent(inout) :: self

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(mpas_field), intent(inout)  :: self
real(kind=kind_real), intent(in) :: zz

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(mpas_field), intent(inout)  :: self
real(kind=kind_real), intent(in) :: zz
type(mpas_field), intent(in)     :: rhs

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

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(mpas_field), intent(inout) :: lhs
type(mpas_field), intent(in)    :: x1
type(mpas_field), intent(in)    :: x2

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(mpas_field), intent(inout) :: fld
type(mpas_field), intent(in)    :: rhs

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
!logical, pointer :: config_do_restart
!integer :: ierr
!real (kind=RKIND), pointer :: dt

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      !clock => fld % domain % clock
!      call mpas_pool_get_config(fld % domain % blocklist % configs, 'config_do_restart', config_do_restart)
!      call mpas_pool_get_config(fld % domain % blocklist % configs, 'config_dt', dt)
      !
      ! If this is a restart run, read the restart stream, else read the input
      ! stream.
      ! Regardless of which stream we read for initial conditions, reset the
      ! input alarms for both input and restart before reading any remaining
      ! input streams.
      !
!      if (config_do_restart) then
!         call MPAS_stream_mgr_read(fld % domain % streamManager, streamID='restart', ierr=ierr)
!      else
!         call MPAS_stream_mgr_read(fld % domain % streamManager, streamID='input', ierr=ierr)
!      end if
!      if (ierr /= MPAS_STREAM_MGR_NOERR) then
!         call mpas_dmpar_global_abort('********************************************************************************', &
!                                      deferredAbort=.true.)
!         call mpas_dmpar_global_abort('Error reading initial conditions', &
!                                      deferredAbort=.true.)
!         call mpas_dmpar_global_abort('********************************************************************************')
!      end if
!      call MPAS_stream_mgr_reset_alarms(fld % domain % streamManager, streamID='input', direction=MPAS_STREAM_INPUT, ierr=ierr)
!      call MPAS_stream_mgr_reset_alarms(fld % domain % streamManager, streamID='restart', direction=MPAS_STREAM_INPUT, ierr=ierr)

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
