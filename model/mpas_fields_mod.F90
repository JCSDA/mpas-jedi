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

use model_oops, only: model_oops_coord
use obsop_apply, only: apply_obsop 
use obsop_parameters, only: compute_parameters
use type_nam, only: namtype
use type_geom, only: geomtype,compute_grid_mesh
use type_odata, only: odatatype
use type_randgen, only: create_randgen
use mpi
use mpas_pool_routines

implicit none
type(mpas_field), intent(in)    :: fld
type(ufo_locs), intent(in)  :: locs
type(mpas_vars), intent(in)     :: vars
type(ufo_geovals), intent(inout) :: gom

logical, save :: interp_initialized = .FALSE.
integer :: n
real(kind=kind_real), parameter :: deg2rad = 0.01745329251
real(kind=kind_real), parameter :: rad2deg = 57.29577954572

integer :: mod_nC,mod_nz,mod_ns
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
real(kind=kind_real), allocatable :: mod_field(:,:)
integer :: obs_num
real(kind=kind_real), allocatable :: obs_lat(:), obs_lon(:) 
real(kind=kind_real), allocatable :: obs_field(:,:)

type(geomtype), save, target :: geom
type(namtype), target :: nam
integer, allocatable :: imask(:)
real(kind=kind_real), allocatable :: area(:),vunit(:)
type(odatatype), save :: odata

integer :: ierr, i, j

type (mpas_pool_type), pointer :: pool_a
type (mpas_pool_iterator_type) :: poolItr
real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b


!Get the Solution dimensions
!---------------------------
mod_nC  = fld%geom%nCells
mod_nz  = fld%geom%nVertLevels
mod_ns  = 1


!Check things are the sizes we expect
!------------------------------------
if( gom%nvar .lt. mod_ns )then
   call abor1_ftn("mpas_fields_mod:interp_tl nvar wrong size")
endif
if( .not. allocated(gom%geovals) )then
   call abor1_ftn("mpas_fields_mod:interp_tl geovals unallocated")
endif
if( size(gom%geovals, 1) .lt. mod_ns )then
   call abor1_ftn("mpas_fields_mod:interp_tl geovals wrong size")
endif
do n=1,mod_ns
   if( allocated(gom%geovals(n)%vals) )then  
      if( gom%geovals(n)%nval .ne. mod_nz )then
         call abor1_ftn("mpas_fields_mod:interp_tl nval wrong size")
      endif
      if( gom%geovals(n)%nobs .ne. mod_nC )then
         call abor1_ftn("mpas_fields_mod:interp_tl nobs wrong size")
      endif
      if( size(gom%geovals(n)%vals, 1) .eq. mod_nz )then
         call abor1_ftn("mpas_fields_mod:interp_tl vals wrong size 1")
      endif
      if( size(gom%geovals(n)%vals, 2) .eq. mod_nC )then
         call abor1_ftn("mpas_fields_mod:interp_tl vals wrong size 2")
      endif       
   endif 
enddo


!Calculate interpolation weight using nicas
!------------------------------------------
if( .NOT. interp_initialized )then
   
   !Initialized flag
   interp_initialized = .TRUE. 

   allocate( mod_lat(mod_nC), mod_lon(mod_nC) )
   mod_lat = fld%geom%latCell
   mod_lon = fld%geom%lonCell

   !*!HARDCODED!*! obsnum to 2 for testing, when deleting change intent to in
   obs_num = 0
!   if (mpp_pe() == 1) then
      obs_num = obs_num + 1
      !Uncomment to give observation on this processor, then put in json file
      !locs%lat(1) = rad2deg*fld%geom%grid_lat(fld%geom%bd%isc+5,fld%geom%bd%jsc+7) + 0.47
      !locs%lon(1) = rad2deg*fld%geom%grid_lon(fld%geom%bd%isc+5,fld%geom%bd%jsc+7) + 0.27
      print*, 'TestOB 1', locs%lat(1), locs%lon(1)
!   elseif (mpp_pe() == 2) then
      obs_num = obs_num + 1
      !Uncomment to give observation on this processor, then put in json file
      !locs%lat(2) = rad2deg*fld%geom%grid_lat(fld%geom%bd%isc+9,fld%geom%bd%jsc+8) + 0.62
      !locs%lon(2) = rad2deg*fld%geom%grid_lon(fld%geom%bd%isc+9,fld%geom%bd%jsc+8) + 0.25
      print*, 'TestOB 2', locs%lat(2), locs%lon(2)
!   endif
   !*!END HARDCODED!*!

   allocate( obs_lat(obs_num), obs_lon(obs_num) )

   !*!HARDCODED!*! test mode
   !Uncomment below lines
!   if (mpp_pe() == 1) then
   obs_lat(1)  = deg2rad * locs%lat(1)
   obs_lon(1)  = deg2rad * locs%lon(1)
!   elseif (mpp_pe() == 2) then
   obs_lat(2)  = deg2rad * locs%lat(2)
   obs_lon(2)  = deg2rad * locs%lon(2)
!   endif
   !*!END HARDCODED!*!

   !obs_lat(:)  = deg2rad * locs%lat(:)
   !obs_lon(:)  = deg2rad * locs%lon(:)

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
   geom%nc0a = mod_nC   ! Number of grid points (local)
   geom%nl0 = 1         ! Number of levels: only one level here (same interpolation for all levels)
   geom%nlev = geom%nl0 ! Copy
   allocate(area(mod_nC))
   allocate(vunit(1))
   allocate(imask(mod_nC))
   area = 1.0           ! Dummy area
   vunit = 1.0          ! Dummy vertical unit
   imask = 1            ! Mask
   call model_oops_coord(geom,mod_lon,mod_lat,area,vunit,imask)

   !Compute grid mesh
   call compute_grid_mesh(nam,geom)

   !Initialize observation operator with observations coordinates (local)
   odata%nobsa = obs_num
   allocate(odata%lonobs(odata%nobsa))
   allocate(odata%latobs(odata%nobsa))
   odata%lonobs = obs_lon
   odata%latobs = obs_lat

   !Setup observation operator
   odata%nam => nam
   odata%geom => geom
   call compute_parameters(odata,.true.)

   deallocate( mod_lat, mod_lon )
   deallocate( obs_lat, obs_lon )

endif


!Make sure the return values are allocated and set
!-------------------------------------------------
do n=1,mod_ns
   if( .not. allocated(gom%geovals(n)%vals) )then
      allocate( gom%geovals(n)%vals(mod_nz,obs_num) )
      gom%geovals(n)%nval = mod_nz
      gom%geovals(n)%nobs = mod_nC
   endif
enddo


!Create Buffer for interpolated values
!--------------------------------------
allocate(mod_field(mod_nC,1))
allocate(obs_field(obs_num,1))

!Interpolate fields to obs locations using pre-calculated weights
!----------------------------------------------------------------
call mpas_pool_begin_iteration(fld % subFields)

do while ( mpas_pool_get_next_member(fld % subFields, poolItr) )
     write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
     ! Pools may in general contain dimensions, namelist options, fields, or other pools,
     ! so we select only those members of the pool that are fields
     if (poolItr % memberType == MPAS_POOL_FIELD) then
     ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
     if (poolItr % dataType == MPAS_POOL_REAL) then
        ! Depending on the dimensionality of the field, we need to set pointers of
        ! the correct type
        if (poolItr % nDims == 1) then
           !call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
           !select case (op_type)
           !case ('TL')
           !   call apply_linop(geov_hinterp_op, r1d_ptr_a(:), fld_dst)
           !case ('AD')
           !   !call apply_linop_ad(hinterp_op,fld_dst,r1d_ptr_a(:))
           !end select
        else if (poolItr % nDims == 2) then
           call mpas_pool_get_array(fld % subFields, trim(poolItr % memberName), r2d_ptr_a)
           !call mpas_pool_get_array(poolItr, trim(poolItr % memberName), r2d_ptr_a)
           do n = 1, mod_nz
              mod_field(:,1) = r2d_ptr_a(n,:)
              call apply_obsop(geom,odata,mod_field,obs_field)
              write(0,*)"n, MIN/MAX: ",minval(r2d_ptr_a(:,n)),maxval(r2d_ptr_a(:,n))
              write(*,*)"n, Interp. value = ", n, obs_field(:,1)
              gom%geovals(1)%vals(n,:) = obs_field(:,1)
           end do
        else if (poolItr % nDims == 3) then
           write(*,*)'Not implemented yet'
           !call abort
        end if
     end if
     end if
end do

!!ua
!do n = 1, mod_nz
!
!   mod_field(:,1) = reshape( fld%Atm%ua(fld%geom%bd%isc:fld%geom%bd%iec,      &
!                                   fld%geom%bd%jsc:fld%geom%bd%jec,      &
!                                   n), [mod_num])
!
!   call apply_obsop(geom,odata,mod_field,obs_field)
!
!   gom%geovals(1)%vals(n,:) = obs_field(:,1)
!enddo
!
!!va
!do n = 1, mod_nz
!
!   mod_field(:,1) = reshape( fld%Atm%va(fld%geom%bd%isc:fld%geom%bd%iec,      &
!                                   fld%geom%bd%jsc:fld%geom%bd%jec,      &
!                                   n), [mod_num])
!
!   call apply_obsop(geom,odata,mod_field,obs_field)
!
!   gom%geovals(2)%vals(n,:) = obs_field(:,1)
!enddo                      
!
!!pt
!do n = 1, mod_nz
!
!   mod_field(:,1) = reshape( fld%Atm%pt(fld%geom%bd%isc:fld%geom%bd%iec,      &
!                                   fld%geom%bd%jsc:fld%geom%bd%jec,      &
!                                   n), [mod_num])
!
!   call apply_obsop(geom,odata,mod_field,obs_field)
!
!   gom%geovals(3)%vals(n,:) = obs_field(:,1)
!enddo 
!
!!q (tracer 1)
!do n = 1, mod_nz
!
!   mod_field(:,1) = reshape( fld%Atm%q(fld%geom%bd%isc:fld%geom%bd%iec,      &
!                                     fld%geom%bd%jsc:fld%geom%bd%jec,      &
!                                     n, 1), [mod_num])
!
!   call apply_obsop(geom,odata,mod_field,obs_field)
!
!   gom%geovals(4)%vals(n,:) = obs_field(:,1)
!enddo 
!
!!delp
!do n = 1, mod_nz
!
!   mod_field(:,1) = reshape( fld%Atm%delp(fld%geom%bd%isc:fld%geom%bd%iec,      &
!                                     fld%geom%bd%jsc:fld%geom%bd%jec,      &
!                                     n), [mod_num])
!
!   call apply_obsop(geom,odata,mod_field,obs_field)
!
!   gom%geovals(5)%vals(n,:) = obs_field(:,1)
!enddo                                                                   

deallocate(mod_field)
deallocate(obs_field)

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
