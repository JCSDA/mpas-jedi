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
use mpas_kinds, only : kind_real
use ioda_locs_mod
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
use mpas_pool_routines
use mpas_constants

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

real(kind=kind_real), parameter :: deg2rad = pii/180.0_kind_real !-BJJ: TODO:  To-be-removed, when MPAS-release updated from Gael.
 
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
    !--- TODO: aux test: BJJ  !- get this from json ???
    nf_aux=19 !15 ! 3
    allocate(fldnames_aux(nf_aux))
    !fldnames_aux(1)='pressure'
    !fldnames_aux(2)='u10'
    !fldnames_aux(3)='t2m'
    !fldnames_aux=(/"pressure", "u10", "t2m"/)
    fldnames_aux = [ character(len=22) :: "pressure", "landmask", "xice", "snowc", "skintemp", "ivgtyp", "isltyp", &
                                          "snowh", "vegfra", "u10", "v10", "lai", "smois", "tslb", "w", &
                                          "index_qc", "index_qi", "re_cloud", "re_ice" ]
                                          !BJJ- "w" is for var_prsi 
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
    call update_mpas_field(self % geom % domain, self % auxFields)
!!    call mpas_pool_get_subpool(self % geom % domain % blocklist % structs,'state',state)
!    call da_fldrms(self % subFields, self % geom % domain % dminfo, prms)
!    allocate(pstat(3, self % nf))
!    call da_gpnorm(self % subFields, self % geom % domain % dminfo, self % nf, pstat)
!    deallocate(pstat)
!!    call abor1_ftn("MPAS test")

     call zeros(self) !-- set zero for self % subFields

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

   temp_filename = config_get_string(c_conf,len(temp_filename),"filename")
   write(*,*)'Reading ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   ! TODO: we can get streamID from json
   !streamID = 'restart'
   !streamID = 'input'
   streamID = 'output'
   !streamID = 'da'
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
   !--TODO: BJJ test. Do I need to "re-calculate"/"update" diagnostic variables ?
   !call update_mpas_field(fld % geom % domain, fld % auxFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % subFields) 
   call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields) 

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)

use duration_mod
   implicit none
   type(mpas_field), intent(inout) :: fld    !< Fields
   type(c_ptr),      intent(in)    :: c_conf !< Configuration
   type(datetime),   intent(inout) :: vdate  !< DateTime
   character(len=20)       :: validitydate
   integer                 :: ierr
   type (MPAS_Time_type)   :: fld_time, write_time
   character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename, temp_filename

   call datetime_to_string(vdate, validitydate)
   write(*,*)'==> write fields at ',trim(validitydate)
   temp_filename = config_get_string(c_conf,len(temp_filename),"filename")
   write(*,*)'==> writing ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(validitydate,dateTimeString,-1)
   ierr = 0
   call mpas_set_time(write_time, dateTimeString=dateTimeString, ierr=ierr)
   fld_time = mpas_get_clock_time(fld % clock, MPAS_NOW, ierr)
   call mpas_get_time(fld_time, dateTimeString=dateTimeString2, ierr=ierr)
   if ( fld_time .NE. write_time ) then
      write(*,*)'write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
      call abor1_ftn('Different times MPAS_stream_mgr_write failed ')
   end if
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)
   fld % manager => fld % geom % domain % streamManager
   ! TODO: we can get streamID from json
   !streamID = 'restart'
   streamID = 'output'
   !streamID = 'da'
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   call da_copy_sub2all_fields(fld % geom % domain, fld % subFields)
   call da_copy_sub2all_fields(fld % geom % domain, fld % auxFields)
   write(*,*)'writing ',trim(filename)
   call mpas_stream_mgr_write(fld % geom % domain % streamManager, streamID=streamID, forceWriteNow=.true., ierr=ierr)
   if ( ierr .ne. 0  ) then
     call abor1_ftn('MPAS_stream_mgr_write failed ierr=',ierr)
   end if

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)

   implicit none
   type(mpas_field),     intent(in)  :: fld
   integer,              intent(in)  :: nf
   real(kind=kind_real), intent(out) :: pstat(3, nf)

   call da_gpnorm(fld % subFields, fld % geom % domain % dminfo, fld%nf, pstat) 

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
   if ((ifdir<1).or.(ifdir>self%nf)) then
      call abor1_ftn("mpas_fields:dirac invalid ifdir")
   endif

   ! Setup Diracs
   call zeros(self)

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
                  r2d_ptr_a(ildir,iCell(idir))=1.0_kind_real
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

   use type_bump, only: bump_type
   use mpas2ufo_vars_mod !, only: usgs_to_crtm_mw, wrf_to_crtm_soil

   implicit none
   type(mpas_field),  intent(in)    :: fld
   type(ioda_locs),   intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom
   
   type(bump_type), pointer :: pbump
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   real(kind=kind_real), allocatable :: tmp_field(:,:)  !< for wspeed/wdir
   
   type (mpas_pool_type), pointer :: pool_ufo  !< pool with ufo variables
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b

   real(kind=kind_real) :: wdir           !< for wind direction
   integer :: ivarw, ivarl, ivari, ivars  !< for sfc fraction indices


   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs 
   write(*,*)'interp: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("nl", fld, locs, vars, gom)
   
   ! Calculate interpolation weight using BUMP
   ! ------------------------------------------
   call initialize_interp(fld%geom, locs, pbump)
   write(*,*)'interp: after initialize_interp'
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
! BJJ: move allocation inside mpas_pool iteration
!   do jvar=1,vars%nv
!      if( .not. allocated(gom%geovals(jvar)%vals) )then
!         gom%geovals(jvar)%nval = fld%geom%nVertLevels
!         allocate( gom%geovals(jvar)%vals(fld%geom%nVertLevels,nobs) )
!         write(*,*) ' gom%geovals(n)%vals allocated'
!      endif
!   enddo
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
   call convert_mpas_field2ufo(fld % geom, fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)

   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then

        if (poolItr % dataType == MPAS_POOL_INTEGER) then
           if (poolItr % nDims == 1) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i1d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = 1
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              mod_field(:,1) = real( i1d_ptr_a(:) )
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(i1d_ptr_a),maxval(i1d_ptr_a)
              call pbump%apply_obsop(mod_field,obs_field)
              gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
           end if
        end if

        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = 1
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              mod_field(:,1) = r1d_ptr_a(:)
              write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r1d_ptr_a),maxval(r1d_ptr_a)
              call pbump%apply_obsop(mod_field,obs_field)
              gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
              write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = fld%geom%nVertLevels
                 if(trim(poolItr % memberName).eq.var_prsi) gom%geovals(ivar)%nval = fld%geom%nVertLevelsP1 !BJJ: Can we do this better ??
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r2d_ptr_a),maxval(r2d_ptr_a)
              do jlev = 1, gom%geovals(ivar)%nval
                 mod_field(:,1) = r2d_ptr_a(jlev,:)
                 call pbump%apply_obsop(mod_field,obs_field)
                 !ORG- gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
                 gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) = obs_field(:,1) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
              end do
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

           else if (poolItr % nDims == 3) then
           end if
        end if

        end if
   end do !- end of pool iteration


   !---add special cases: var_sfc_wspeed and/or var_sfc_wdir
   if ( (ufo_vars_getindex(vars,var_sfc_wspeed)    .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_wdir) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_wspeed and/or var_sfc_wdir'

     !- allocate
     allocate(tmp_field(nobs,2))

     !- read/interp.
     call mpas_pool_get_array(fld % auxFields, "u10", r1d_ptr_a)
     write(*,*) 'MIN/MAX of u10=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(:)
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,1)=obs_field(:,1)
     call mpas_pool_get_array(fld % auxFields, "v10", r1d_ptr_a)
     write(*,*) 'MIN/MAX of v10=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(:)
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,2)=obs_field(:,1)

     !- allocate geoval & put values for var_sfc_wspeed
     ivar = ufo_vars_getindex(vars, var_sfc_wspeed)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       gom%geovals(ivar)%vals(1,:) = sqrt( tmp_field(:,1)**2 + tmp_field(:,2)**2 ) ! ws = sqrt(u**2+v**2) [m/s]
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wspeed),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_wdir
     ivar = ufo_vars_getindex(vars, var_sfc_wdir)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       do ii=1,nobs
         call uv_to_wdir(tmp_field(ii,1), tmp_field(ii,2), wdir) ! uu, vv, wind10_direction in radian
         gom%geovals(ivar)%vals(1,ii) = wdir / deg2rad           ! radian -> degree
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wdir),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- deallocate
     deallocate(tmp_field)
    
   endif  !---end special cases


   !---add special cases: var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp
   if ( (ufo_vars_getindex(vars,var_sfc_landtyp)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_vegtyp)  .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_soiltyp) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_landtyp, var_sfc_vegtyp, or var_sfc_soiltyp'

     call mpas_pool_get_array(fld % auxFields, "ivgtyp", i1d_ptr_a)
     write(*,*) 'MIN/MAX of ivgtyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)

     call mpas_pool_get_array(fld % auxFields, "isltyp", i1d_ptr_b)
     write(*,*) 'MIN/MAX of isltyp=',minval(i1d_ptr_b),maxval(i1d_ptr_b)

     !- allocate geoval & put values for var_sfc_landtyp
     ivar = ufo_vars_getindex(vars, var_sfc_landtyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       !BJJ: How I implement the nearest neighbor interp. using existing barycentric interp. weights
       !write(*,*) 'BJJ iobs, n_s, n_src, n_dst =', ii, pbump%obsop%h%n_s, pbump%obsop%h%n_src, pbump%obsop%h%n_dst
       !write(*,*) '          size(row), size(col) =', size(pbump%obsop%h%row), size(pbump%obsop%h%col)
       !do ii=1,odata%h%n_s;write(*,*) 'BJJ ith operator, ROW = ROW + S * COL =', ii, odata%h%row(ii), odata%h%S(ii), odata%h%col(ii)
       !compare odata%h%S( 3*(ii-1) + 1, 3*(ii-1) + 2, 3*(ii-1) + 3 ), find minimum, then specify odata%h%col( that index ) for geoval.
       !write(*,*) minloc(odata%h%S), minloc(odata%h%S( 3*(ii-1)+1:3*(ii-1)+3 ))
       do ii=1,nobs
         !write(*,*) ' i-th obs, max weight, Cell index=',ii, maxval(odata%h%S( 3*(ii-1)+1:3*(ii-1)+3 )), &
         !            maxloc(odata%h%S( 3*(ii-1)+1:3*(ii-1)+3 )), 3*(ii-1) + maxloc(odata%h%S( 3*(ii-1)+1:3*(ii-1)+3 )), &
         !            odata%h%col( 3*(ii-1) + maxloc(odata%h%S( 3*(ii-1)+1:3*(ii-1)+3 )) )
         jj=3*(ii-1) + maxloc(pbump%obsop%h%S( 3*(ii-1)+1:3*(ii-1)+3 ),1) !nearest-interp. / maximum-weight specified.
         gom%geovals(ivar)%vals(1,ii) = i1d_ptr_a(jj)
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_landtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_vegtyp
     ivar = ufo_vars_getindex(vars, var_sfc_vegtyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       do ii=1,nobs
         jj=3*(ii-1) + maxloc(pbump%obsop%h%S( 3*(ii-1)+1:3*(ii-1)+3 ),1)
         gom%geovals(ivar)%vals(1,ii) = convert_type_veg( i1d_ptr_a(jj) )  !nearest-interp. / maximum-weight specific.
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_vegtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_soiltyp
     ivar = ufo_vars_getindex(vars, var_sfc_soiltyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       do ii=1,nobs
         jj=3*(ii-1) + maxloc(pbump%obsop%h%S( 3*(ii-1)+1:3*(ii-1)+3 ),1)
         gom%geovals(ivar)%vals(1,ii) = convert_type_soil( i1d_ptr_b(jj) )  !nearest-interp. / maximum-weight specific.
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_soiltyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

   endif  !---end special cases


   !---add special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac
   !---    simple interpolation now, but can be more complex: Consider FOV ??
   if ( (ufo_vars_getindex(vars,var_sfc_wfrac)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_lfrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_ifrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_sfrac) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, or var_sfc_sfrac'

     call mpas_pool_get_array(fld % auxFields, "landmask", i1d_ptr_a) !"land-ocean mask (1=land ; 0=ocean)"
     write(*,*) 'MIN/MAX of landmask=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
     call mpas_pool_get_array(fld % auxFields, "xice", r1d_ptr_a)     !"fractional area coverage of sea-ice"
     write(*,*) 'MIN/MAX of xice=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
     call mpas_pool_get_array(fld % auxFields, "snowc", r1d_ptr_b)    !"flag for snow on ground (=0 no snow; =1,otherwise"
     write(*,*) 'MIN/MAX of snowc=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
     write(*,*) 'MIN/MAX of lnad+xice+snowc=', &
                minval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b),maxval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b)

     ivarw = ufo_vars_getindex(vars, var_sfc_wfrac)
     ivarl = ufo_vars_getindex(vars, var_sfc_lfrac)
     ivari = ufo_vars_getindex(vars, var_sfc_ifrac)
     ivars = ufo_vars_getindex(vars, var_sfc_sfrac)

     !--- Land first. will be adjusted later
     ivar = ivarl
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = real(i1d_ptr_a(:))
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- determine ICE
     ivar = ivari
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = r1d_ptr_a(:)
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_ifrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- detemine/adjust SNOW & SEA
     ivar = ivars
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       gom%geovals(ivarw)%nval = 1
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       allocate( gom%geovals(ivarw)%vals(gom%geovals(ivarw)%nval,nobs) )
       mod_field(:,1) = r1d_ptr_b(:)
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       do ii=1, nobs
         if(gom%geovals(ivari)%vals(1,ii).gt.0.0_kind_real) then
           gom%geovals(ivar)%vals(1,ii) = min( gom%geovals(ivar)%vals(1,ii), 1.0_kind_real - gom%geovals(ivari)%vals(1,ii) )
           gom%geovals(ivarw)%vals(1,ii)= 1.0_kind_real - gom%geovals(ivari)%vals(1,ii) - gom%geovals(ivar)%vals(1,ii)
         else
           gom%geovals(ivarw)%vals(1,ii)= 1.0_kind_real - gom%geovals(ivarl)%vals(1,ii)
         endif
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_sfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wfrac),minval(gom%geovals(ivarw)%vals),maxval(gom%geovals(ivarw)%vals)
     endif
     !--- Final adjust LAND
     ivar = ivarl
     if(ivar .ne. -1) then
       do ii=1, nobs
         gom%geovals(ivar)%vals(1,ii) = max( 1.0_kind_real - gom%geovals(ivarw)%vals(1,ii) - gom%geovals(ivari)%vals(1,ii) &
                                                          - gom%geovals(ivars)%vals(1,ii), 0.0_kind_real)
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !do ii=17,19  !1,nobs
     !write(*,*) gom%geovals(ivarl)%vals(1,ii), gom%geovals(ivarw)%vals(1,ii), &
     !           gom%geovals(ivari)%vals(1,ii), gom%geovals(ivars)%vals(1,ii)
     !enddo
   endif  !---end special cases

   !--- OMG: adjust between sfc coverage & sfc type
   !         see wrfda/da_get_innov_vector_crtm.inc#L521
   if ( (ufo_vars_getindex(vars,var_sfc_wfrac)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_lfrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_ifrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_sfrac) .ne. -1) ) then
   if ( (ufo_vars_getindex(vars,var_sfc_landtyp)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_vegtyp)  .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_soiltyp) .ne. -1) ) then
     do ii=1, nobs
       if(gom%geovals(ivarl)%vals(1,ii) .gt. 0.0_kind_real) then
         if(nint(gom%geovals(ufo_vars_getindex(vars,var_sfc_soiltyp))%vals(1,ii)) .eq. 9 .or. &
            nint(gom%geovals(ufo_vars_getindex(vars,var_sfc_vegtyp))%vals(1,ii)) .eq. 13 ) then
           gom%geovals(ivari)%vals(1,ii) = min( gom%geovals(ivari)%vals(1,ii) + gom%geovals(ivarl)%vals(1,ii), 1.0_kind_real )
           gom%geovals(ivarl)%vals(1,ii) = 0.0_kind_real
         endif
       endif
     enddo
   endif
   endif  !--- OMG: end

   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

   write(*,*) '---- Leaving interp ---'
end subroutine interp

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, vars, gom)

   use type_bump, only: bump_type
   
   implicit none
   !type(mpas_field),  intent(in)    :: fld
   type(mpas_field),  intent(inout) :: fld
   type(ioda_locs),   intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom
   
   type(bump_type), pointer :: pbump
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   
   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   !--- BJJ tmp: Trajectory hack
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   type (mpas_pool_type), pointer :: subFields1, auxFields1
   type (mpas_pool_type), pointer :: Traj_subFields, Traj_auxFields
   type (mpas_pool_type), pointer :: subFields3, auxFields3
   integer :: ierr=0
   
   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs 
   write(*,*)'interp_tl: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("tl", fld, locs, vars, gom)
   
   ! Calculate interpolation weight using BUMP
   ! ------------------------------------------
   call initialize_interp(fld%geom, locs, pbump)
   write(*,*)'interp_tl: after initialize_interp'
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
   do jvar=1,vars%nv
      if( .not. allocated(gom%geovals(jvar)%vals) )then
         gom%geovals(jvar)%nval = fld%geom%nVertLevels
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
   write(0,*)'interp_tl: vars%nv       : ',vars%nv
   write(0,*)'interp_tl: vars%fldnames : ',vars%fldnames
   
!--- BJJ tmp: Trajectory hack
!---------Save Original TL--
   call mpas_pool_create_pool(subFields1)
   call mpas_pool_clone_pool(fld % subFields, subFields1)
   call mpas_pool_create_pool(auxFields1)
   call mpas_pool_clone_pool(fld % auxFields, auxFields1)
!---------Read Traj. w/ fld--
   streamID = 'output'
   fld % manager => fld % geom % domain % streamManager
   filename='/home/vagrant/build/mpas-bundle/mpas/test/x1.2562.init.nc'
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(*,*)'Reading ',trim(filename)
   call MPAS_stream_mgr_read(fld % manager, streamID=streamID, &
                           & rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      call abor1_ftn('MPAS_stream_mgr_read failed ierr=',ierr)
   end if
   !-- BJJ test. Do I need to "re-calculate"/"update" diagnostic variables ?
   !call update_mpas_field(fld % geom % domain, fld % auxFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % subFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields)
!--------Save Traj 
   call mpas_pool_create_pool(Traj_subFields)
   call mpas_pool_clone_pool(fld % subFields, Traj_subFields)
   call mpas_pool_create_pool(Traj_auxFields)
   call mpas_pool_clone_pool(fld % auxFields, Traj_auxFields)
!--------Recover Original TL : ONLY NEED subFields, NOT auxFields
   call mpas_pool_empty_pool(fld % subFields)
   call mpas_pool_destroy_pool(fld % subFields)
   fld % nf = 5
   call mpas_pool_create_pool(fld % subFields,fld % nf)
   call mpas_pool_clone_pool(subFields1, fld % subFields)
!   call mpas_pool_clone_pool(auxFields1, fld % auxFields)
!--------Remove temporary pool
   call mpas_pool_empty_pool(subFields1)
   call mpas_pool_destroy_pool(subFields1)
   call mpas_pool_empty_pool(auxFields1)
   call mpas_pool_destroy_pool(auxFields1)
!--------------------------
   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   !call convert_mpas_field2ufoTL(fld % subFields, fld % auxFields, pool_b, vars % fldnames, vars % nv) !--pool_b is new pool with ufo_vars
   call convert_mpas_field2ufoTL(Traj_subFields, Traj_auxFields, &
        fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              write(*,*) "interp_tl: var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex: ivar=",ivar
              do jlev = 1, gom%geovals(ivar)%nval
                 mod_field(:,1) = r2d_ptr_a(jlev,:)
                 call pbump%apply_obsop(mod_field,obs_field)
                 !ORG- gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
                 gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) = obs_field(:,1) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do


   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

!--------Remove temporary pool
   call mpas_pool_empty_pool(Traj_subFields)
   call mpas_pool_destroy_pool(Traj_subFields)

   write(*,*) '---- Leaving interp_tl ---'
end subroutine interp_tl

! ------------------------------------------------------------------------------

subroutine interp_ad(fld, locs, vars, gom)

   use type_bump, only: bump_type

   implicit none
   type(mpas_field),  intent(inout) :: fld
   type(ioda_locs),   intent(in)    :: locs
   type(ufo_vars),    intent(in)    :: vars
   type(ufo_geovals), intent(inout) :: gom

   type(bump_type), pointer :: pbump

   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)

   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   !--- BJJ tmp: Trajectory hack
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   type (mpas_pool_type), pointer :: subFields1, auxFields1
   type (mpas_pool_type), pointer :: Traj_subFields, Traj_auxFields
   type (mpas_pool_type), pointer :: subFields3, auxFields3
   integer :: ierr=0

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCells
   nobs = locs%nlocs

   call interp_checks("ad", fld, locs, vars, gom)

   ! Calculate interpolation weight using BUMP
   ! ------------------------------------------
   call initialize_interp(fld%geom, locs, pbump)
   write(*,*)'interp_ad: after initialize_interp'

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))

   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'interp_ad: vars%nv       : ',vars%nv
   write(0,*)'interp_ad: vars%fldnames : ',vars%fldnames

!--- BJJ tmp: Trajectory hack
!---------Save Original AD--
   call mpas_pool_create_pool(subFields1)
   call mpas_pool_clone_pool(fld % subFields, subFields1)
   call mpas_pool_create_pool(auxFields1)
   call mpas_pool_clone_pool(fld % auxFields, auxFields1)
!---------Read Traj. w/ fld--
   streamID = 'output'
   fld % manager => fld % geom % domain % streamManager
   filename='/home/vagrant/build/mpas-bundle/mpas/test/x1.2562.init.nc'
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(*,*)'Reading ',trim(filename)
   call MPAS_stream_mgr_read(fld % manager, streamID=streamID, &
                           & rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      call abor1_ftn('MPAS_stream_mgr_read failed ierr=',ierr)
   end if
   !-- BJJ test. Do I need to "re-calculate"/"update" diagnostic variables ?
   !call update_mpas_field(fld % geom % domain, fld % auxFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % subFields)
   call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields)
!--------Save Traj 
   call mpas_pool_create_pool(Traj_subFields)
   call mpas_pool_clone_pool(fld % subFields, Traj_subFields)
   call mpas_pool_create_pool(Traj_auxFields)
   call mpas_pool_clone_pool(fld % auxFields, Traj_auxFields)
!--------Recover Original AD : ONLY NEED subFields, NOT auxFields
    write(*,*) 'here ---Recover Original AD'
   call mpas_pool_empty_pool(fld % subFields)
   call mpas_pool_destroy_pool(fld % subFields)
   fld % nf = 5
   call mpas_pool_create_pool(fld % subFields,fld % nf)
   call mpas_pool_clone_pool(subFields1, fld % subFields)
!   call mpas_pool_clone_pool(auxFields1, fld % auxFields)
!--------Remove temporary pool
   call mpas_pool_empty_pool(subFields1)
   call mpas_pool_destroy_pool(subFields1)
   call mpas_pool_empty_pool(auxFields1)
   call mpas_pool_destroy_pool(auxFields1)
!--------------------------
   !BJJ tmp:call convert_mpas_field2ufoTL(fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars
   call convert_mpas_field2ufoTL(Traj_subFields, Traj_auxFields, &
        fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              r2d_ptr_a=0.0
              write(*,*) "Interp. var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex, ivar=",ivar
              do jlev = 1, gom%geovals(ivar)%nval
                 !ORG- obs_field(:,1) = gom%geovals(ivar)%vals(jlev,:)
                 obs_field(:,1) = gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
                 call pbump%apply_obsop_ad(obs_field,mod_field)
                 r2d_ptr_a(jlev,:) = 0.0_kind_real
                 r2d_ptr_a(jlev,:) = r2d_ptr_a(jlev,:) + mod_field(:,1)
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do

   !call convert_mpas_field2ufoAD(fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars
   call convert_mpas_field2ufoAD(Traj_subFields, Traj_auxFields, &
        fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

!--------Remove temporary pool
   call mpas_pool_empty_pool(Traj_subFields)
   call mpas_pool_destroy_pool(Traj_subFields)

   write(*,*) '---- Leaving interp_ad ---' 
end subroutine interp_ad

! ------------------------------------------------------------------------------

subroutine initialize_interp(grid, locs, pbump)

   use mpas_geom_mod, only: mpas_geom
   use type_bump, only: bump_type
   
   implicit none
   type(mpas_geom),          intent(in)  :: grid
   type(ioda_locs),          intent(in)  :: locs
   type(bump_type), pointer, intent(out) :: pbump
   
   logical, save :: interp_initialized = .FALSE.
   type(bump_type), save, target :: bump
   
   integer :: mod_nz,mod_num,obs_num
   real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
   
   real(kind=kind_real), allocatable :: area(:),vunit(:,:)
   logical, allocatable :: lmask(:,:)
   
   integer :: ii, jj, ji, jvar, jlev
  
   !Get the Solution dimensions
   !---------------------------
   mod_nz  = grid%nVertLevels
   mod_num = grid%nCells
   obs_num = locs%nlocs 
   write(*,*)'initialize_interp mod_num,obs_num = ',mod_num,obs_num
   
   !Calculate interpolation weight using BUMP
   !------------------------------------------
   if (.NOT.interp_initialized) then
      allocate( mod_lat(mod_num), mod_lon(mod_num) )
      mod_lat = grid%latCell / deg2rad !- to Degrees
      mod_lon = grid%lonCell / deg2rad !- to Degrees
   
      !Important namelist options
      bump%nam%prefix       = 'oops_data'  ! Prefix for files output
      bump%nam%nobs         = obs_num      ! Number of observations
      bump%nam%obsop_interp = 'bilin'      ! Interpolation type (bilinear)
      bump%nam%obsdis       = 'local'      ! Observation distribution parameter ('random', 'local', or 'adjusted')
   
      !Less important namelist options (should not be changed)
      bump%nam%default_seed        = .true.
      bump%nam%new_hdiag           = .false.
      bump%nam%new_param           = .false.
      bump%nam%check_adjoints      = .false.
      bump%nam%check_pos_def       = .false.
      bump%nam%check_sqrt          = .false.
      bump%nam%check_dirac         = .false.
      bump%nam%check_randomization = .false.
      bump%nam%check_consistency   = .false.
      bump%nam%check_optimality    = .false.
      bump%nam%new_lct             = .false.
      bump%nam%new_obsop           = .true.

      !Initialize geometry
      allocate(area(mod_num))
      allocate(vunit(mod_num,1))
      allocate(lmask(mod_num,1))
      area  = 1.0          ! Dummy area, unit [m^2]
      vunit = 1.0          ! Dummy vertical unit
      lmask = .true.       ! Mask

      !Initialize BUMP
      call bump%setup_online(mpi_comm_world,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                             nobs=obs_num,lonobs=locs%lon(:),latobs=locs%lat(:))

      !Release memory
      deallocate(area)
      deallocate(vunit)
      deallocate(lmask)
      deallocate( mod_lat, mod_lon )

      interp_initialized = .TRUE. 
   endif

   pbump => bump

end subroutine initialize_interp

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)

   implicit none
   character(len=2),  intent(in) :: cop
   type(mpas_field),  intent(in) :: fld
   type(ioda_locs),   intent(in) :: locs
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
   if (cop/="tl" .and. cop/='nl') then !BJJ or cop='ad'.   why only check ad??? because ad is set/allocated in UFO side ??
                                                           ! also, the dimension is not general for all possible geovals.
   if (.not.gom%linit) then
      call abor1_ftn(cinfo//"geovals not initialized")
   endif
   do jvar=1,vars%nv
      if (allocated(gom%geovals(jvar)%vals)) then  
         if( gom%geovals(jvar)%nval .ne. fld%geom%nVertLevels )then
            write(*,*) jvar, gom%geovals(jvar)%nval, fld%geom%nVertLevels
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
   
   integer :: nmga,jC,jl,nf
   integer,allocatable :: imask(:,:)
   real(kind=kind_real),allocatable :: lon(:),lat(:),area(:),vunit(:,:)
   real(kind=kind_real) :: sigmaup,sigmadn
   
   type (mpas_pool_iterator_type) :: poolItr
   !real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   !real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   !real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: idx_var
   
   ! Define local number of gridpoints
   nmga = self%geom%nCells
   
   ! Allocation
   allocate(lon(nmga))
   allocate(lat(nmga))
   allocate(area(nmga))
   allocate(vunit(nmga,self%geom%nVertLevels))
   allocate(imask(nmga,self%geom%nVertLevels))
   
   ! Copy coordinates
   do jC=1,self%geom%nCells
     lon(jC) = self%geom%lonCell(jC) / deg2rad !- to Degrees
     lat(jC) = self%geom%latCell(jC) / deg2rad !- to Degrees
     area(jC) = self%geom%areaCell(jC)
   enddo

   imask = 1
   
   ! Define vertical unit
   do jl=1,self%geom%nVertLevels
     !sigmaup = self%geom%ak(jl+1)/101300.0+self%geom%bk(jl+1) ! si are now sigmas
     !sigmadn = self%geom%ak(jl  )/101300.0+self%geom%bk(jl  )
     !vunit(jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
     vunit(:,jl) = real(jl,kind=kind_real)
   enddo

   ! Should this come from self/vars?
   nf = self%nf
   write(0,*)'convert_to_ug: nf=',nf
   !if (.not. self%geom%hydrostatic) nf = 7

   ! Create unstructured grid
   call create_unstructured_grid(ug, nmga, self%geom%nVertLevels, nf, 1, lon, lat, area, vunit, imask)

   !! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
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
   
   integer :: jC,jl

   type (mpas_pool_iterator_type) :: poolItr
   !real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   !real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   !real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: idx_var

   ! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
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
