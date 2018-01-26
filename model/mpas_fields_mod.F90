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
use mpi ! only MPI_COMM_WORLD
use mpas_stream_manager
use type_mpl, only: mpl,mpl_recv,mpl_send,mpl_bcast

implicit none
private

public :: mpas_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp_tl, interp_ad, define_ug, convert_to_ug, convert_from_ug
public :: mpas_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold MPAS fields
!type :: mpas_field
!  type (domain_type), pointer :: domain 
!  type (core_type), pointer :: corelist
!  !type (dm_info), pointer :: dminfo
!  type (mpas_geom), pointer :: geom                                 !< Number of unstructured grid cells
!  integer :: nf                           !< Number of variables in fld
!  !integer :: sum_scalar                   !< Number of variables in fld
!  !integer :: sum_aero                     !< Number of variables in fld
!  !integer :: ns                           !< Number of surface fields (x1d [nCells])
!  character(len=22), allocatable  :: fldnames(:)      !< Variable identifiers
!  type (mpas_pool_type), pointer  :: subFields
!end type mpas_field

#include "mpas_fields_oops_type.inc"

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

    use mpas_kind_types

    implicit none

    type(mpas_field), intent(inout) :: self
    !type(mpas_geom),  intent(inout), pointer    :: geom
    type(mpas_geom),  intent(in), pointer    :: geom
    type(mpas_vars),  intent(in)    :: vars

    !type (dm_info), pointer :: dminfo
    integer :: nsize
    integer :: mpi_comm
    real(kind=kind_real), allocatable :: pstat(:, :)
    !type (dm_info) :: dminfo

    ! from the namelist
    self % nf = 1 ! vars % nv
    allocate(self % fldnames(self % nf))
    self % fldnames(1) = "theta_m" !vars % fldnames(:)
    write(*,*)'allocate ',self % fldnames(:)
   
    ! coming from mpas_subdriver
    write(*,*)'Before calling the subdriver of mpas'
    !if ( geom % use_mpi ) then
    !if ( self % domain % dminfo % initialized_mpi  ) then
    !write(*,*)'dminfo % initialized_mpi: ',dminfo % initialized_mpi,geom % use_mpi
    !if ( geom % use_mpi  ) then
       !self % domain % dminfo => dminfo
       call mpas_init(self % corelist, self % domain, mpi_comm=MPI_COMM_WORLD)
    !else
    !   call mpas_init(self % corelist, self % domain)
    !   geom % use_mpi = .true.
    !   dminfo % initialized_mpi = .true.
    !end if
    ! dminfo % initialized_mpi
    !write(*,*)'After calling the subdriver of mpas'

    ! update geom
    if (associated(geom)) then
      write(*,*)'geom associated'
      self % geom => geom
    else
      write(*,*)'geom not associated'
      self % geom =>null()
    end if
    ! This or create a subpool dminfo and clone
!    geom % dminfo => self % domain % dminfo

    ! Create a subpool from allfields
!    nsize = da_common_vars(self % domain % blocklist % allFields, self % fldnames)
!    if ( self % nf .ne. nsize  ) then
!       call abor1_ftn("mpas_fields:create: dimension mismatch ",self % nf, nsize)
!    end  if
    !write(0,*)'-- Create a sub Pool from list of variable ',nsize
!    call mpas_pool_create_pool(self % subFields,self % nf)
    call da_make_subpool(self % domain % blocklist % allFields, self % subFields, self % fldnames)
    allocate(pstat(3, self % nf))
    call da_gpnorm(self % subFields, self % nf, pstat)
    write(0,*)'Pstat: ',pstat
    deallocate(pstat)
    !call abor1_ftn('lfric_fields_mod: create: call not implemented yet')
    write(*,*)'Leaving create'

    
 
    return

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(mpas_field), intent(inout) :: self
type (MPAS_streamManager_type), pointer :: manager
integer :: ierr = 0 
  
   call mpas_timer_set_context( self % domain )
   
   !
   ! Set up run-time streams
   !
   !   call MPAS_stream_mgr_init(self % domain % streamManager, self % domain % ioContext, & 
   !      self % domain % clock, self % domain % blocklist % allFields, &
   !      self % domain % packages, self % domain % blocklist % allStructs)
   !
   !   call add_stream_attributes(domain)
   !
   !   ierr = domain % core % setup_immutable_streams(domain % streamManager)
   !   if ( ierr /= 0 ) then
   !      call mpas_dmpar_global_abort('ERROR: Immutable streams setup failed for core ' // trim(domain % core % coreName))
   !   end if
   !
   !   mgr_p = c_loc(domain_ptr % streamManager)
   !   call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
   !   if (c_ierr /= 0) then
   !      call mpas_dmpar_abort(domain_ptr % dminfo)
   !   end if

   if (.not. associated(manager)) then  
     write(0,*)'GD manager not associated'
   else
     write(0,*)'GD manager associated'
   end if
   !  manager % ioContext => self % domain % ioContext
   !  manager % allFields => self % domain % blocklist % allFields
   !  manager % allPackages => self % domain % Packages
   !  manager % allStructs => self % domain % blocklist % allStructs
   !  manager % streamClock => self % domain % clock
   !  manager % numStreams = 0
   !  manager % errorLevel = MPAS_STREAM_ERR_SILENT

 
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   write(*,*)'--> deallocate subFields Pool'
   call mpas_pool_empty_pool(self % subFields)
   call mpas_pool_destroy_pool(self % subFields)
   write(*,*)'--> deallocate domain and core'
   call mpas_finalize(self % corelist, self % domain)
   write(*,*)'--> done'

   return

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
real(kind=kind_real) :: zz
   
   zz = 10.
   write(0,*)'Calling da_setval: ' 
   call da_setval(self % subFields, zz)
   !call da_random(self % subFields)

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
   
   ! Duplicate the members of rhs into self and do a deep copy
   ! of the fields from self % subFields to rhs % subFields
   call mpas_pool_empty_pool(self % subFields)
   call mpas_pool_destroy_pool(self % subFields)
   self % nf = rhs % nf
   call mpas_pool_create_pool(self % subFields,self % nf)
   !write(0,*)'copy pool before 1' 
   call mpas_pool_clone_pool(rhs % subFields, self % subFields)
   !call mpas_pool_clone_pool(rhs % domain % blocklist % allFields, self % domain % blocklist % allFields)
   !write(0,*)'copy pool before 2' 
   ! We should consider adding a subroutine just updating the fields
   ! call mpas_pool_copy_fied() 

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

return

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

implicit none
type(mpas_field), intent(inout) :: fld      !< Fields
type(c_ptr), intent(in)          :: c_conf   !< Configuration
type(datetime), intent(inout)    :: vdate    !< DateTime
character (len=StrKIND) :: dateTimeString
character(len=20) :: sdate

write(*,*)'inside read'
sdate = config_get_string(c_conf,len(sdate),"date")
! GD look at oops/src/util/datetime_mod.F90
! we probably need to extract from vdate a string to enforce the reading ..
! and then can be like this ....
!dateTimeString = '2010-10-24_04:00:00'
!write(0,*)''
!write(0,*)'Reading ',dateTimeString
!call MPAS_stream_mgr_read(field0 % domain % streamManager,streamID='restart',when=dateTimeString)
sdate = config_get_string(c_conf,len(sdate),"date")
call datetime_set(sdate, vdate)
write(*,*)'outside read'

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
real(kind=kind_real), intent(out) :: pstat(3, nf)
real(kind=kind_real), allocatable :: pstat2(:, :)

allocate(pstat2(3, fld% nf))
pstat(3, nf) = 15.0
write(0,*)'Inside gpnorm 1'
call da_gpnorm(fld % subFields, fld%nf, pstat2) 
write(0,*)'Inside gpnorm 2'
deallocate(pstat2)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
implicit none
type(mpas_field), intent(in)      :: fld
real(kind=kind_real), intent(out) :: prms
real(kind=kind_real), allocatable :: pstat(:,:)

allocate(pstat(3,fld % nf))
pstat(:, :) = 15.0
write(0,*)'Inside fldrms 1'
call da_gpnorm(fld % subFields, fld % nf, pstat) 
write(0,*)'Inside fldrms 2'
prms =  pstat(2,1)
deallocate(pstat)

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

subroutine define_ug(self, ug)
use unstructured_grid_mod
implicit none
type(mpas_field), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug

integer,allocatable :: imask(:,:)
real(kind=kind_real),allocatable :: lon(:),lat(:),area(:),vunit(:)
!integer :: nc0a,ic0a,jx,jy,jl,jf,joff
integer :: nc0a,ic0a,jC,jl,jf,joff

! Define local index
nc0a = 0
do jC=1,self%geom%nCells
  if (self%geom%iproc(jC)==mpl%myproc) nc0a = nc0a+1
enddo

! Allocation
allocate(lon(nc0a))
allocate(lat(nc0a))
allocate(area(nc0a))
allocate(vunit(self%geom%nVertLevels))
allocate(imask(nc0a,self%geom%nVertLevels))

! Copy coordinates
ic0a = 0
do jC=1,self%geom%nCells
  if (self%geom%iproc(jC)==mpl%myproc) then
    ic0a = ic0a+1
    lon(ic0a) = self%geom%lonCell(jC)
    lat(ic0a) = self%geom%latCell(jC)
    area(ic0a) = self%geom%areaCell(jC)
  endif
enddo
imask = 1

! Define vertical unit
do jl=1,self%geom%nVertLevels
  vunit(jl) = real(jl,kind=kind_real)
enddo

! Create unstructured grid
call create_unstructured_grid(ug, nc0a, self%geom%nVertLevels, self%nf, 1, lon, lat, area, vunit, imask)

end subroutine define_ug

! ------------------------------------------------------------------------------

subroutine convert_to_ug(self, ug)
use unstructured_grid_mod
implicit none
type(mpas_field), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug

integer :: ic0a,jC,jl,jf,joff

! Copy field
ic0a = 0
do jC=1,self%geom%nCells
  if (self%geom%iproc(jC)==mpl%myproc) then
    ic0a = ic0a+1
    do jf=1,self%nf
      joff = (jf-1)*self%geom%nVertLevels
      do jl=1,self%geom%nVertLevels
!        W/ mpas field type
!        ug%fld(ic0a,jl,jf,1) = self%gfld3d(jx,jy,joff+jl)
      enddo
    enddo
  endif
enddo

end subroutine convert_to_ug

! ------------------------------------------------------------------------------

subroutine convert_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(mpas_field), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug

integer :: ic0a,jC,jl,jf,joff,nbuf,jbuf,iproc
real(kind=kind_real),allocatable :: rbuf(:),sbuf(:)

! Copy field
ic0a = 0
do jC=1,self%geom%nCells
  if (self%geom%iproc(jC)==mpl%myproc) then
    ic0a = ic0a+1
    do jf=1,self%nf
      joff = (jf-1)*self%geom%nVertLevels
      do jl=1,self%geom%nVertLevels
!        W/ mpas field type
!        self%gfld3d(jx,jy,joff+jl) = ug%fld(ic0a,jl,jf,1)
      enddo
    enddo
  endif
enddo

! Communication
if (mpl%main) then
   do iproc=1,mpl%nproc
      if (iproc/=mpl%ioproc) then
         ! Allocation
         nbuf = count(self%geom%iproc==iproc)*self%geom%nVertLevels*self%nf
         allocate(rbuf(nbuf))

         ! Receive data on ioproc
         call mpl_recv(nbuf,rbuf,iproc,mpl%tag)

         ! Format data
         jbuf = 0
         do jC=1,self%geom%nCells
           if (self%geom%iproc(jC)==iproc) then
             do jf=1,self%nf
               joff = (jf-1)*self%geom%nVertLevels
               do jl=1,self%geom%nVertLevels
                 jbuf = jbuf+1
!                 W/ mpas field type
!                 self%gfld3d(jx,jy,joff+jl) = rbuf(jbuf)
               enddo
             enddo
           endif
         enddo

         ! Release memory
         deallocate(rbuf)
      end if
   end do
else
   ! Allocation
   nbuf = count(self%geom%iproc==mpl%myproc)*self%nf*self%geom%nVertLevels
   allocate(sbuf(nbuf))

   ! Format data
   jbuf = 0
   do jC=1,self%geom%nCells
     if (self%geom%iproc(jC)==mpl%myproc) then
       do jf=1,self%nf
         joff = (jf-1)*self%geom%nVertLevels
         do jl=1,self%geom%nVertLevels
           jbuf = jbuf+1
!           W/ mpas field type
!           sbuf(jbuf) = self%gfld3d(jx,jy,joff+jl)
         enddo
       enddo
     endif
   enddo

   ! Send data to ioproc
   call mpl_send(nbuf,sbuf,mpl%ioproc,mpl%tag)

   ! Release memory
   deallocate(sbuf)
end if
mpl%tag = mpl%tag+1

! Broadcast
!W/ mpas field type
!call mpl_bcast(self%gfld3d,mpl%ioproc)

end subroutine convert_from_ug

! ------------------------------------------------------------------------------

end module mpas_fields_mod
