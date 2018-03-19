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
use ufo_vars_mod
use mpas_kinds
use ufo_locs_mod
use ufo_geovals_mod

use mpas_derived_types
use mpas_framework
use mpas_kind_types
!use init_atm_core_interface
use mpas_subdriver
use atm_core
use mpas4da_mod
use mpas2ufo_vars_mod
use mpi ! only MPI_COMM_WORLD
use mpas_stream_manager
use type_mpl, only: mpl,mpl_recv,mpl_send,mpl_bcast
use mpas_pool_routines

implicit none
private

public :: mpas_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp, interp_tl, interp_ad, convert_to_ug, convert_from_ug, &
        & dirac
public :: mpas_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold MPAS fields
type :: mpas_field
  !type (domain_type), pointer :: domain          ! NOW in geom, For convenience not used 
  !type (core_type), pointer :: corelist          ! NOW in geom, For convenience not used
  type (mpas_geom), pointer :: geom              ! grid and MPI infos
  integer :: nf                                  ! Number of variables in fld
  character(len=22), allocatable  :: fldnames(:) ! Variable identifiers
  type (mpas_pool_type), pointer  :: subFields   !---> state variables (to be analyzed)
  type (mpas_pool_type), pointer  :: auxFields   !---> auxiliary variables, such as pressure, t2m, u10, v10, Tsfc
  type (MPAS_streamManager_type), pointer :: manager
  type (MPAS_Clock_type), pointer :: clock
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

    use mpas_kind_types

    implicit none

    type(mpas_field), intent(inout)       :: self
    type(mpas_geom),  intent(in), pointer :: geom
    type(ufo_vars),   intent(in)          :: vars

    integer :: nsize, nfields
    integer :: ierr

    !aux test: BJJ
    integer :: nf_aux
    character(len=22), allocatable  :: fldnames_aux(:)

    !-- fortran level test (temporally sit here)
    real(kind=kind_real), allocatable :: pstat(:, :)
    real(kind=kind_real)              :: prms
    type (MPAS_Time_type) :: local_time, write_time, fld_time
    character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename
    character (len=StrKIND) :: dateTimeString_oops


    ! from the namelist
    self % nf =  vars % nv
    allocate(self % fldnames(self % nf))
    self % fldnames(:) = vars % fldnames(:)
    write(*,*)'self % nf =',self % nf
    write(*,*)'allocate ::',self % fldnames(:)
   
    ! link geom
    if (associated(geom)) then
      self % geom => geom
    else
      write(*,*)'mpas_fields: geom not associated'
      call abor1_ftn("mpas_fields: geom not associated")
    end if

    write(0,*)'-- Create a sub Pool from list of variable ',self % nf
    call da_make_subpool(self % geom % domain, self % subFields, self % nf, self % fldnames, nfields)
    if ( self % nf .ne. nfields  ) then
       call abor1_ftn("mpas_fields:create: dimension mismatch ",self % nf, nfields)
    end  if
    !--- aux test: BJJ
    nf_aux=3
    allocate(fldnames_aux(nf_aux))
    fldnames_aux(1)='pressure'
    fldnames_aux(2)='u10'
    fldnames_aux(3)='t2m'
    write(0,*)'-- Create a sub Pool for auxFields'
    call da_make_subpool(self % geom % domain, self % auxFields, nf_aux, fldnames_aux, nfields)
    if ( nf_aux .ne. nfields  ) then
       call abor1_ftn("mpas_fields:create: dimension mismatch ",self % nf, nfields)
    end  if

    ! clock creation
    allocate(self % clock)
    call atm_simulation_clock_init(self % clock, self % geom % domain % blocklist % configs, ierr)
    if ( ierr .ne. 0 ) then
       call abor1_ftn("mpas_fields: atm_simulation_clock_init problem")
    end if

    !-------------------------------------------------------------
    ! Few temporary tests
    !-------------------------------------------------------------
    call update_mpas_field(self % geom % domain, self % subFields)
!    call mpas_pool_get_subpool(self % geom % domain % blocklist % structs,'state',state)
    call da_fldrms(self % subFields, self % geom % domain % dminfo, prms)
    allocate(pstat(3, self % nf))
    call da_gpnorm(self % subFields, self % geom % domain % dminfo, self % nf, pstat)
    deallocate(pstat)
!    call abor1_ftn("MPAS test")

   ! mpas_get_time(curr_time, YYYY, MM, DD, DoY, H, M, S, S_n, S_d, dateTimeString, ierr)
   ! call mpas_get_time(curr_time, YYYY, MM, DD, H, M, S, dateTimeString, ierr)
   ! call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
   ! call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)
   ! currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
   ! startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
   ! currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
   ! xtimeTime = currTime - startTime
   ! call mpas_get_timeInterval(interval=xtimeTime, S=s, S_n=s_n, S_d=s_d, ierr=ierr)
   ! xtime_s = (s + s_n / s_d)


!--------------------------------------------------------------------------------------- 

!   write(0,*)''
!   write(0,*)'=========== MPAS WRITE A RESTART ===================='

   !call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
!   call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
!   all mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
!   call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
!   call atm_compute_restart_diagnostics(state, 1, diag, mesh)
 
    return

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

   implicit none
   type(mpas_field), intent(inout) :: self
   integer :: ierr = 0 
  
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   write(*,*)'--> deallocate subFields Pool'
   if (associated(self % subFields)) then
      call mpas_pool_empty_pool(self % subFields)
      call mpas_pool_destroy_pool(self % subFields)
   end if
   call mpas_destroy_clock(self % clock, ierr)
   if ( ierr .ne. 0  ) then
      write(*,*)'mpas_fields eallocate clock failed'
   end if
   write(*,*)'--> mpas_fields done deallocate'

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
   !call da_setval(self % subFields, zz)
   call da_random(self % subFields)

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
   call mpas_pool_clone_pool(rhs % subFields, self % subFields)
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
   type(mpas_field),     intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz

   call da_self_mult(self % subFields, zz)

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

   implicit none
   type(mpas_field),     intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz
   type(mpas_field),     intent(in)    :: rhs

   call da_axpy(self % subFields, rhs % subFields, zz)

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)

   implicit none
   type(mpas_field),     intent(in)    :: fld1, fld2
   real(kind=kind_real), intent(inout) :: zprod

   call da_dot_product(fld1 % subFields, fld2 % subFields, fld1 % geom % domain % dminfo, zprod)

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
   type(c_ptr),      intent(in)    :: c_conf   !< Configuration
   type(datetime),   intent(inout) :: vdate    !< DateTime
   character(len=20)       :: sdate
   type (MPAS_Time_type)   :: local_time
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   integer                 :: ierr = 0

   write(*,*)'==> read fields'
   sdate = config_get_string(c_conf,len(sdate),"date")
   call datetime_set(sdate, vdate)
#define READ_TEST_READY
#ifdef READ_TEST_READY
   temp_filename = config_get_string(c_conf,len(temp_filename),"filename")
   write(*,*)'Reading ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   !streamID = 'restart'
   !streamID = 'input'
   streamID = 'output'
   ierr = 0
   fld % manager => fld % geom % domain % streamManager
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(sdate,dateTimeString,1)
   write(*,*)'dateTimeString: ',trim(dateTimeString)
   call mpas_set_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
   !call mpas_set_clock_time(fld % clock, local_time, MPAS_START_TIME)
   call mpas_set_clock_time(fld % geom % domain % clock, local_time, MPAS_START_TIME)
   call mpas_expand_string(dateTimeString, -1, temp_filename, filename)
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(*,*)'Reading ',trim(filename)
   call MPAS_stream_mgr_read(fld % manager, streamID=streamID, &
                           & when=dateTimeString, rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      call abor1_ftn('MPAS_stream_mgr_read failed ierr=',ierr)
   end if
   !-- BJJ test. Do I need to "re-calculate"/"update" diagnostic variables ?
   !call update_mpas_field(fld % geom % domain, fld % subFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % subFields) 
   call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields) 
#endif

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)

   implicit none
   type(mpas_field), intent(inout) :: fld    !< Fields
   type(c_ptr),      intent(in)    :: c_conf !< Configuration
   type(datetime),   intent(inout) :: vdate  !< DateTime
   character(len=20)       :: sdate
   integer                 :: ierr
   type (MPAS_Time_type)   :: fld_time, write_time
   character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename, temp_filename

#define WRITE_TEST_READY
#ifdef WRITE_TEST_READY
   sdate = config_get_string(c_conf,len(sdate),"date")
   write(*,*)'==> write fields ',trim(sdate)
   call datetime_set(sdate,vdate)
   temp_filename = config_get_string(c_conf,len(temp_filename),"filename")
   write(*,*)'Writing ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....

   !dateTimeString_oops = '2010-10-23T00.00.00Z'  ! <- sdate form
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(sdate,dateTimeString,-1)
   !write(*,*) 'dateTimeString_oops=',dateTimeString_oops
   !write(*,*) 'dateTimeString=',dateTimeString
   ierr = 0
   call mpas_set_time(write_time, dateTimeString=dateTimeString, ierr=ierr)
   fld_time = mpas_get_clock_time(fld % clock, MPAS_START_TIME, ierr)
   call mpas_get_time(fld_time, dateTimeString=dateTimeString2, ierr=ierr)
   if ( fld_time .NE. write_time ) then
      write(*,*)'write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
      call abor1_ftn('Different times MPAS_stream_mgr_write failed ')
   end if
   !temp_filename = 'restart.out.$Y-$M-$D_$h.$m.$s.nc'
   write(*,*) 'temp_filename = ',trim(temp_filename)
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)
   write(*,*) '     filename = ',trim(filename)
   fld % manager => fld % geom % domain % streamManager
   !streamID = 'restart'
   streamID = 'output'
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   call da_copy_sub2all_fields(fld % geom % domain, fld % subFields)
   call da_copy_sub2all_fields(fld % geom % domain, fld % auxFields)
   write(*,*)'writing ',trim(filename)
   call mpas_stream_mgr_write(fld % geom % domain % streamManager, streamID=streamID, forceWriteNow=.true., ierr=ierr)
   if ( ierr .ne. 0  ) then
     call abor1_ftn('MPAS_stream_mgr_write failed ierr=',ierr)
   end if
#endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)

   implicit none
   type(mpas_field),     intent(in)  :: fld
   integer,              intent(in)  :: nf
   real(kind=kind_real), intent(out) :: pstat(3, nf)

   write(0,*)'Inside gpnorm 1'
   call da_gpnorm(fld % subFields, fld % geom % domain % dminfo, fld%nf, pstat) 
   write(0,*)'Inside gpnorm 2'

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)

   implicit none
   type(mpas_field),     intent(in)  :: fld
   real(kind=kind_real), intent(out) :: prms

   call da_fldrms(fld % subFields, fld % geom % domain % dminfo, prms)

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)

use iso_c_binding
use mpas_pool_routines

   implicit none
   type(mpas_field), intent(inout) :: self
   type(c_ptr),      intent(in)    :: c_conf   !< Configuration
   integer                :: ndir,idir,ildir,ifdir,itiledir
   integer,allocatable    :: iCell(:)
   character(len=3)       :: idirchar
   character(len=StrKIND) :: dirvar
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a

   ! Get Diracs positions
   ndir = config_get_int(c_conf,"ndir")
   allocate(iCell(ndir))

   do idir=1,ndir
      write(idirchar,'(i3)') idir
      iCell(idir) = config_get_int(c_conf,"iCell("//trim(adjustl(idirchar))//")")
   end do
   ildir = config_get_int(c_conf,"ildir")
   ifdir = config_get_int(c_conf,"ifdir")
   dirvar = config_get_string(c_conf,len(dirvar),"dirvar")

   ! Check
   if (ndir<1) call abor1_ftn("mpas_fields:dirac non-positive ndir")
   if (any(iCell<1).or.any(iCell>self%geom%nCells)) then
      call abor1_ftn("mpas_fields:dirac invalid iCell")
   endif
   if ((ildir<1).or.(ildir>self%geom%nVertLevels)) then
      call abor1_ftn("mpas_fields:dirac invalid ildir")
   endif
   if ((ifdir<1).or.(ifdir>5)) then
      call abor1_ftn("mpas_fields:dirac invalid ifdir")
   endif

   ! Setup Diracs
   call zeros(self)
   !do idir=1,ndir
   !   self%fld(iCell(idir),ildir,ifdir) = 1.0
   !end do

   call mpas_pool_begin_iteration(self % subFields)

   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
          else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
              if( trim(dirvar) .eq. trim(poolItr % memberName) ) then
                do idir=1, ndir
                  r2d_ptr_a(ildir,iCell(idir))=1.0
                end do
                write(*,*) ' Dirac is set in ',ndir,'locations for',trim(poolItr % memberName)
              end if
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine interp(fld, locs, vars, gom)

   use obsop_apply, only: apply_obsop 
   use type_geom, only: geomtype
   use type_odata, only: odatatype
   
   implicit none
   type(mpas_field),  intent(in)    :: fld
   type(ufo_locs),    intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom
   
   type(geomtype), pointer :: pgeom
   type(odatatype), pointer :: odata
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   
   type (mpas_pool_type), pointer :: pool_b
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   
   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs 
   write(*,*)'interp: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("tl", fld, locs, vars, gom)
   
   ! Calculate interpolation weight using nicas
   ! ------------------------------------------
   call initialize_interp(fld%geom, locs, pgeom, odata)
   write(*,*)'interp: after initialize_interp'
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
   do jvar=1,vars%nv
      if( .not. allocated(gom%geovals(jvar)%vals) )then
         gom%geovals(jvar)%nval = fld%geom%nVertLevels
         gom%geovals(jvar)%nobs = ngrid
         allocate( gom%geovals(jvar)%vals(fld%geom%nVertLevels,nobs) )
         write(*,*) ' gom%geovals(n)%vals allocated'
      endif
   enddo
   gom%linit = .true.
   

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))
   
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'interp: vars%nv       : ',vars%nv
   write(0,*)'interp: vars%fldnames : ',vars%fldnames
   

   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufo(fld % subFields, fld % auxFields, pool_b, vars % fldnames, vars % nv) !--pool_b is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_b)
   do while ( mpas_pool_get_next_member(pool_b, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_a)
              write(*,*) "interp: var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex: ivar=",ivar
              do jlev = 1, fld%geom%nVertLevels
                 mod_field(:,1) = r2d_ptr_a(jlev,:)
                 call apply_obsop(pgeom,odata,mod_field,obs_field)
                 gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
              end do

           else if (poolItr % nDims == 3) then
           end if
        end if
        end if
   end do


   deallocate(mod_field)
   deallocate(obs_field)
   
end subroutine interp

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, vars, gom)

   use obsop_apply, only: apply_obsop 
   use type_geom, only: geomtype
   use type_odata, only: odatatype
   
   implicit none
   type(mpas_field),  intent(in)    :: fld
   type(ufo_locs),    intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom
   
   type(geomtype), pointer :: pgeom
   type(odatatype), pointer :: odata
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   
   type (mpas_pool_type), pointer :: pool_b
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   
   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs 
   write(*,*)'interp_tl: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("tl", fld, locs, vars, gom)
   
   ! Calculate interpolation weight using nicas
   ! ------------------------------------------
   call initialize_interp(fld%geom, locs, pgeom, odata)
   write(*,*)'interp_tl: after initialize_interp'
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
   do jvar=1,vars%nv
      if( .not. allocated(gom%geovals(jvar)%vals) )then
         gom%geovals(jvar)%nval = fld%geom%nVertLevels
         gom%geovals(jvar)%nobs = ngrid
         allocate( gom%geovals(jvar)%vals(fld%geom%nVertLevels,nobs) )
         write(*,*) ' gom%geovals(n)%vals allocated'
      endif
   enddo
   gom%linit = .true.
   

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))
   
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'interp: vars%nv       : ',vars%nv
   write(0,*)'interp: vars%fldnames : ',vars%fldnames
   

   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufoTL(fld % subFields, fld % auxFields, pool_b, vars % fldnames, vars % nv) !--pool_b is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_b)
   do while ( mpas_pool_get_next_member(pool_b, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_a)
              write(*,*) "interp_tl: var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex: ivar=",ivar
              do jlev = 1, fld%geom%nVertLevels
                 mod_field(:,1) = r2d_ptr_a(jlev,:)
                 call apply_obsop(pgeom,odata,mod_field,obs_field)
                 gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do


   deallocate(mod_field)
   deallocate(obs_field)
   
end subroutine interp_tl

! ------------------------------------------------------------------------------

subroutine interp_ad(fld, locs, vars, gom)

   use obsop_apply, only: apply_obsop_ad
   use type_geom, only: geomtype
   use type_odata, only: odatatype

   implicit none
   type(mpas_field),  intent(inout) :: fld
   type(ufo_locs),    intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom

   type(geomtype), pointer :: pgeom
   type(odatatype), pointer :: odata

   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)

   type (mpas_pool_type), pointer :: pool_b
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs
   write(*,*)'interp_ad: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("ad", fld, locs, vars, gom)

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))

   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'interp_ad: vars%nv       : ',vars%nv
   write(0,*)'interp_ad: vars%fldnames : ',vars%fldnames

   call convert_mpas_field2ufoTL(fld % subFields, fld % auxFields, pool_b, vars % fldnames, vars % nv) !--pool_b is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_b)
   do while ( mpas_pool_get_next_member(pool_b, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_a)
              r2d_ptr_a=0.0
              write(*,*) "Interp. var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex, ivar=",ivar
              do jlev = 1, fld%geom%nVertLevels
!TL                 mod_field(:,1) = r2d_ptr_a(jlev,:)
!TL                 call apply_obsop(pgeom,odata,mod_field,obs_field)
!TL                 gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
                 obs_field(:,1) = gom%geovals(ivar)%vals(jlev,:)
                 call apply_obsop_ad(pgeom,odata,obs_field,mod_field)
                 r2d_ptr_a(jlev,:) = r2d_ptr_a(jlev,:) + mod_field(:,1)
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do

   call convert_mpas_field2ufoAD(fld % subFields, fld % auxFields, pool_b, vars % fldnames, vars % nv) !--pool_b is new pool with ufo_vars

   deallocate(mod_field)
   deallocate(obs_field)

   
end subroutine interp_ad

! ------------------------------------------------------------------------------

subroutine initialize_interp(grid, locs, pgeom, pdata)

   use mpas_geom_mod, only: mpas_geom
   use model_oops, only: model_oops_coord
   use obsop_parameters, only: compute_parameters
   use type_nam, only: namtype
   use type_geom, only: geomtype,compute_grid_mesh
   use type_odata, only: odatatype
   use type_randgen, only: create_randgen
   
   implicit none
   type(mpas_geom),          intent(in)  :: grid
   type(ufo_locs),           intent(in)  :: locs
   type(geomtype), pointer,  intent(out) :: pgeom
   type(odatatype), pointer, intent(out) :: pdata
   
   logical, save :: interp_initialized = .FALSE.
   type(geomtype), save, target :: geom
   type(odatatype), save, target :: odata
   
   integer :: mod_nz,mod_num,obs_num
   real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
   
   type(namtype), target :: nam
   integer, allocatable :: imask(:)
   real(kind=kind_real), allocatable :: area(:),vunit(:)
   
   integer :: ii, jj, ji, jvar, jlev
   
   !Get the Solution dimensions
   !---------------------------
   mod_nz  = grid%nVertLevels
   mod_num = grid%nCells
   obs_num = locs%nlocs 
   write(*,*)'initialize_interp mod_num,obs_num = ',mod_num,obs_num
   
   !Calculate interpolation weight using nicas
   !------------------------------------------
   if (.NOT.interp_initialized) then
      allocate( mod_lat(mod_num), mod_lon(mod_num) )
      mod_lat = grid%latCell
      mod_lon = grid%lonCell
   
      !Important namelist options
      nam%obsop_interp = 'bilin' ! Interpolation type (bilinear)
      nam%obsdis = 0.0           ! Observation distribution parameter (0.0 => local distribution, 1.0 => perfect load balancing)
      nam%datadir = '.'          ! Data directory
      nam%prefix = 'oops_data'   ! Prefix for files output
   
      !Less important namelist options (should not be changed)
      nam%default_seed = .true.
      nam%model = 'oops'
      nam%mask_type = 'none'
      nam%new_hdiag = .false.
      nam%displ_diag = .false.
      nam%new_param = .false.
      nam%new_lct = .false.
      nam%mask_check = .false.
      nam%new_obsop = .true.
      nam%check_dirac = .false.
      nam%nc3 = 1
      nam%dc = 1.0
   
      !Initialize random number generator
      call create_randgen(nam)
   
      !Initialize geometry
      geom%nc0a = mod_num  ! Number of grid points (local)
      geom%nl0 = 1         ! Number of levels: only one level here (same interpolation for all levels)
      geom%nlev = geom%nl0 ! Copy
      allocate(area(mod_num))
      allocate(vunit(1))
      allocate(imask(mod_num))
      area = 1.0           ! Dummy area
      vunit = 1.0          ! Dummy vertical unit
      imask = 1            ! Mask
      call model_oops_coord(geom,mod_lon,mod_lat,area,vunit,imask)
      deallocate(area)
      deallocate(vunit)
      deallocate(imask)
   
      !Compute grid mesh
      call compute_grid_mesh(nam,geom)
   
      !Initialize observation operator with observations coordinates (local)
      odata%nobsa = obs_num
      allocate(odata%lonobs(odata%nobsa))
      allocate(odata%latobs(odata%nobsa))
      odata%lonobs(:) = locs%lon(:)
      odata%latobs(:) = locs%lat(:)
   
      !Setup observation operator
      odata%nam => nam
      odata%geom => geom
      call compute_parameters(odata,.true.)
   
      deallocate( mod_lat, mod_lon )

      interp_initialized = .TRUE. 
   endif

   pgeom => geom
   pdata => odata

end subroutine initialize_interp

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)

   implicit none
   character(len=2),  intent(in) :: cop
   type(mpas_field),  intent(in) :: fld
   type(ufo_locs),    intent(in) :: locs
   type(ufo_vars),    intent(in) :: vars
   type(ufo_geovals), intent(in) :: gom
   
   integer :: jvar
   character(len=26) :: cinfo
   
   cinfo="mpas_fields:checks "//cop//" : "
   
   !Check things are the sizes we expect
   !------------------------------------
   if (gom%nobs /= locs%nlocs ) then
      call abor1_ftn(cinfo//"geovals wrong size")
   endif
   if( gom%nvar .ne. vars%nv )then
      call abor1_ftn(cinfo//"nvar wrong size")
   endif
   if( .not. allocated(gom%geovals) )then
      call abor1_ftn(cinfo//"geovals unallocated")
   endif
   if( size(gom%geovals) .ne. vars%nv )then
      call abor1_ftn(cinfo//"geovals wrong size")
   endif
   if (cop/="tl") then
   if (.not.gom%linit) then
      call abor1_ftn(cinfo//"geovals not initialized")
   endif
   do jvar=1,vars%nv
      if (allocated(gom%geovals(jvar)%vals)) then  
         if( gom%geovals(jvar)%nval .ne. fld%geom%nVertLevels )then
            call abor1_ftn(cinfo//"nval wrong size")
         endif
         if( gom%geovals(jvar)%nobs .ne. locs%nlocs )then
            call abor1_ftn(cinfo//"nobs wrong size")
         endif
         if( size(gom%geovals(jvar)%vals, 1) .ne. fld%geom%nVertLevels )then
            call abor1_ftn(cinfo//"vals wrong size 1")
         endif
         if( size(gom%geovals(jvar)%vals, 2) .ne. locs%nlocs )then
            call abor1_ftn(cinfo//"vals wrong size 2")
         endif       
      else
        call abor1_ftn(cinfo//"vals not allocated")
      endif 
   enddo
   endif
   
   write(*,*)'interp_checks ',cinfo,' done'
   
end subroutine interp_checks

! ------------------------------------------------------------------------------

subroutine convert_to_ug(self, ug)

   use mpas_pool_routines
   use unstructured_grid_mod
   
   implicit none
   type(mpas_field),        intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   
   integer :: nc0a,ic0a,jC,jl,nf
   integer,allocatable :: imask(:,:)
   real(kind=kind_real),allocatable :: lon(:),lat(:),area(:),vunit(:)
   real(kind=kind_real) :: sigmaup,sigmadn
   
   type (mpas_pool_iterator_type) :: poolItr
   !real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   !real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   !real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: idx_var
   
   ! Define local number of gridpoints
   nc0a = self%geom%nCells
   
   ! Allocation
   allocate(lon(nc0a))
   allocate(lat(nc0a))
   allocate(area(nc0a))
   allocate(vunit(self%geom%nVertLevels))
   allocate(imask(nc0a,self%geom%nVertLevels))
   
   ! Copy coordinates
   do jC=1,self%geom%nCells
     lon(jC) = self%geom%lonCell(jC)
     lat(jC) = self%geom%latCell(jC)
     area(jC) = self%geom%areaCell(jC)
   enddo
write(0,*)'BJJ area MIN/MAX: ',minval(area),maxval(area)

   imask = 1
   
   ! Define vertical unit
   do jl=1,self%geom%nVertLevels
     !sigmaup = self%geom%ak(jl+1)/101300.0+self%geom%bk(jl+1) ! si are now sigmas
     !sigmadn = self%geom%ak(jl  )/101300.0+self%geom%bk(jl  )
     !vunit(jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
     vunit(jl) = real(jl,kind=kind_real)
   enddo

   ! Should this come from self/vars?
   nf = self%nf !5
   write(0,*)'convert_to_ug: nf=',nf
   !if (.not. self%geom%hydrostatic) nf = 7

   ! Create unstructured grid
   call create_unstructured_grid(ug, nc0a, self%geom%nVertLevels, nf, 1, lon, lat, area, vunit, imask)

   !! Copy field
   !do jC=1,self%geom%nCells
   !  do jl=1,self%geom%nVertLevels
   !!     ug%fld(jC,jl,1,1) = self%Atm%ua(jx,jy,jl)
   !!     ug%fld(jC,jl,2,1) = self%Atm%va(jx,jy,jl)
   !!     ug%fld(jC,jl,3,1) = self%Atm%pt(jx,jy,jl)
   !!     ug%fld(jC,jl,4,1) = self%Atm%q(jx,jy,jl,1)
   !!     ug%fld(jC,jl,5,1) = self%Atm%delp(jx,jy,jl)
   !!     if (.not. self%geom%hydrostatic) ug%fld(jC,jl,6,1) = self%Atm%w(jx,jy,jl)
   !!     if (.not. self%geom%hydrostatic) ug%fld(jC,jl,7,1) = self%Atm%delz(jx,jy,jl)
   !  enddo
   !enddo

   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
   
              idx_var=-999
              if(trim(poolItr % memberName).eq.'theta'   ) idx_var=1
              if(trim(poolItr % memberName).eq.'rho'     ) idx_var=2
              if(trim(poolItr % memberName).eq.'index_qv') idx_var=3
              if(trim(poolItr % memberName).eq.'uReconstructZonal') idx_var=4
              if(trim(poolItr % memberName).eq.'uReconstructMeridional') idx_var=5
              if(idx_var.gt.0) then
                 write(*,*) '  sub. convert_to_ug, poolItr % memberName=',trim(poolItr % memberName)
                 do jC=1,self%geom%nCells
                   do jl=1,self%geom%nVertLevels
                     ug%fld(jC,jl,idx_var,1) = r2d_ptr_a(jl,jC)
                   enddo
                 enddo
              endif
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

end subroutine convert_to_ug

! -----------------------------------------------------------------------------

subroutine convert_from_ug(self, ug)

   use mpas_pool_routines
   use unstructured_grid_mod

   implicit none
   type(mpas_field),        intent(inout) :: self
   type(unstructured_grid), intent(in)    :: ug
   
   integer :: ic0a,jC,jl

   type (mpas_pool_iterator_type) :: poolItr
   !real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   !real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   !real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: idx_var

   ! Copy field
   !ic0a = 0
   !do jy=self%geom%bd%jsc,self%geom%bd%jec
   !  do jx=self%geom%bd%isc,self%geom%bd%iec
   !    ic0a = ic0a+1
   !    do jl=1,self%geom%nlevs
   !        self%Atm%ua(jx,jy,jl) = ug%fld(ic0a,jl,1,1)
   !        self%Atm%va(jx,jy,jl) = ug%fld(ic0a,jl,2,1)
   !        self%Atm%pt(jx,jy,jl) = ug%fld(ic0a,jl,3,1)
   !        self%Atm%q(jx,jy,jl,1) = ug%fld(ic0a,jl,4,1)
   !        self%Atm%delp(jx,jy,jl) = ug%fld(ic0a,jl,5,1)
   !        if (.not. self%geom%hydrostatic) self%Atm%w(jx,jy,jl) = ug%fld(ic0a,jl,6,1)
   !        if (.not. self%geom%hydrostatic) self%Atm%delz(jx,jy,jl) = ug%fld(ic0a,jl,7,1)
   !    enddo
   !  enddo
   !enddo

   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)

              idx_var=-999
              if(trim(poolItr % memberName).eq.'theta'   ) idx_var=1
              if(trim(poolItr % memberName).eq.'rho'     ) idx_var=2
              if(trim(poolItr % memberName).eq.'index_qv') idx_var=3
              if(trim(poolItr % memberName).eq.'uReconstructZonal') idx_var=4
              if(trim(poolItr % memberName).eq.'uReconstructMeridional') idx_var=5
              if(idx_var.gt.0) then
                 write(*,*) '  sub. convert_from_ug, poolItr % memberName=',trim(poolItr % memberName)
                 do jC=1,self%geom%nCells
                   do jl=1,self%geom%nVertLevels
                     r2d_ptr_a(jl,jC) = ug%fld(jC,jl,idx_var,1)
                   enddo
                 enddo
              end if
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

end subroutine convert_from_ug

! ------------------------------------------------------------------------------

end module mpas_fields_mod
