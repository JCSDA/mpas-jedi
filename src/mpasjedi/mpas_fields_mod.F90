! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpas_fields_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use string_utils, only: swap_name_member
use unstructured_interpolation_mod

! saber
use interpolatorbump_mod,         only: bump_interpolator

!ufo
use ufo_vars_mod, only: MAXVARLEN, ufo_vars_getindex
use ufo_locations_mod
use ufo_geovals_mod, only: ufo_geovals

!MPAS-Model
use atm_core, only: atm_simulation_clock_init, atm_compute_output_diagnostics
use mpas_constants
use mpas_derived_types
use mpas_kind_types, only: StrKIND
use mpas_pool_routines
use mpas_stream_manager
use mpas_timekeeping

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas4da_mod
use mpas2ufo_vars_mod, only: w_to_q, theta_to_temp


private

public :: mpas_fields, mpas_fields_registry, &
          interp_checks, &
          create_fields, delete_fields, &
          copy_fields, copy_pool, &
          update_diagnostic_fields, &
          mpas_hydrometeor_fields,  &
          mpas_re_fields, &
          cellCenteredWindFields, &
          moistureFields, &
          analysisThermoFields, &
          modelThermoFields

! ------------------------------------------------------------------------------

   !> Fortran derived type to hold MPAS field
   type :: mpas_fields
     private

     type (mpas_geom), pointer, public :: geom                            ! grid and MPI infos
     type (MPAS_streamManager_type), pointer, public :: manager
     type (MPAS_Clock_type), pointer, public :: clock
     integer, public :: nf                                                ! Number of variables in subFields
     character(len=MAXVARLEN), allocatable, public :: fldnames(:)         ! Variable identifiers
     type (mpas_pool_type), pointer, public        :: subFields => null() !---> state variables (to be analyzed)
     integer, public :: nf_ci                                             ! Number of variables in CI
     character(len=MAXVARLEN), allocatable, public :: fldnames_ci(:)      ! Control increment identifiers

     contains

     procedure :: axpy         => axpy_
     procedure :: dot_prod     => dot_prod_
     procedure :: gpnorm       => gpnorm_
     procedure :: random       => random_
     procedure :: rms          => rms_
     procedure :: self_add     => self_add_
     procedure :: self_schur   => self_schur_
     procedure :: self_mult    => self_mult_
     procedure :: self_sub     => self_sub_
     procedure :: zeros        => zeros_
     procedure :: ones         => ones_

     procedure :: change_resol => change_resol_fields
     procedure :: copy         => copy_fields
     procedure :: create       => create_fields
     procedure :: delete       => delete_fields
     procedure :: read_file    => read_fields
     procedure :: write_file   => write_fields
     procedure :: serial_size  => serial_size
     procedure :: serialize    => serialize_fields
     procedure :: deserialize  => deserialize_fields

     generic :: has => has_field, has_fields
     procedure :: has_field
     procedure :: has_fields

   end type mpas_fields

!   abstract interface
!
!   ! ------------------------------------------------------------------------------
!
!      subroutine read_file_(self, f_conf, vdate)
!         import mpas_fields, fckit_configuration, datetime
!         implicit none
!         class(mpas_fields),         intent(inout) :: self
!         type(fckit_configuration), intent(in)    :: f_conf
!         type(datetime),            intent(inout) :: vdate
!      end subroutine read_file_
!
!   ! ------------------------------------------------------------------------------
!
!   end interface

   character(len=MAXVARLEN) :: mpas_hydrometeor_fields(6) = &
      [ character(len=MAXVARLEN) :: &
      "index_qc", "index_qi", "index_qr", &
      "index_qs", "index_qg", "index_qh" ]
   character(len=MAXVARLEN) :: mpas_re_fields(3) = &
      [ character(len=MAXVARLEN) :: &
      "re_cloud", "re_ice  ", "re_snow " ]
   character(len=MAXVARLEN), parameter :: cellCenteredWindFields(2) = &
      [character(len=MAXVARLEN) :: &
       'uReconstructZonal', 'uReconstructMeridional']
   character(len=MAXVARLEN), parameter :: moistureFields(2) = &
      [character(len=MAXVARLEN) :: &
       'index_qv', 'spechum']
   character(len=MAXVARLEN), parameter :: analysisThermoFields(2) = &
      [character(len=MAXVARLEN) :: &
       'surface_pressure', 'temperature']
   character(len=MAXVARLEN), parameter :: modelThermoFields(4) = &
      [character(len=MAXVARLEN) :: &
       'index_qv', 'pressure', 'rho', 'theta']


   integer, parameter    :: max_string=8000
   character(max_string) :: message

#define LISTED_TYPE mpas_fields

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_fields_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include <oops/util/linkedList_c.f>

! ------------------------------------------------------------------------------

subroutine create_fields(self, geom, vars, vars_ci)

    implicit none

    class(mpas_fields),   intent(inout)       :: self
    type(mpas_geom),      intent(in), pointer :: geom
    type(oops_variables), intent(in)          :: vars, vars_ci

    integer :: ivar, ierr

    self % nf = vars % nvars()
    allocate(self % fldnames(self % nf))
    do ivar = 1, self % nf
       self % fldnames(ivar) = trim(vars % variable(ivar))
    end do

    self % nf_ci = vars_ci % nvars()
    allocate(self % fldnames_ci(self % nf_ci))
    do ivar = 1, self % nf_ci
       self % fldnames_ci(ivar) = trim(vars_ci % variable(ivar))
    end do

    write(message,*) "DEBUG: create_fields: self % fldnames(:) =",self % fldnames(:)
    call fckit_log%debug(message)

    ! link geom
    if (associated(geom)) then
      self % geom => geom
    else
      call abor1_ftn("--> create_fields: geom not associated")
    end if

    ! clock creation
    allocate(self % clock)
    call atm_simulation_clock_init(self % clock, self % geom % domain % blocklist % configs, ierr)
    if ( ierr .ne. 0 ) then
       call abor1_ftn("--> create_fields: atm_simulation_clock_init problem")
    end if

!    write(*,*)'--> create_fields: sub Pool from list of variable ',self % nf
    call create_pool(self % geom % domain, self % nf, self % fldnames, self % subFields)

    call self%zeros() !-- set zero for self % subFields

    return

end subroutine create_fields

! ------------------------------------------------------------------------------

subroutine create_pool(domain, nf_in, fldnames, pool)

    implicit none
    type (domain_type), pointer, intent(in) :: domain
    integer, intent(in) :: nf_in
    character (len=*), intent(in) :: fldnames(:)
    type (mpas_pool_type), pointer, intent(out) :: pool

    integer :: nfields

    call da_make_subpool(domain, pool, nf_in, fldnames, nfields)

    if ( nf_in .ne. nfields  ) then
       write(message,*) "--> create_pool: dimension mismatch ", nf_in, nfields
       call abor1_ftn(message)
    end  if

end subroutine create_pool

! ------------------------------------------------------------------------------

subroutine delete_fields(self)

   implicit none
   class(mpas_fields), intent(inout) :: self
   integer :: ierr = 0

   if (allocated(self % fldnames)) deallocate(self % fldnames)
   if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)

!   write(*,*)'--> delete_fields: deallocate subFields Pool'
   call delete_pool(self % subFields)

   call mpas_destroy_clock(self % clock, ierr)
   if ( ierr .ne. 0  ) then
      write(*,*) '--> delete_fields deallocate clock failed'
   end if
!   write(*,*)'--> delete_fields done'

   return

end subroutine delete_fields

! ------------------------------------------------------------------------------

subroutine delete_pool(pool)

   implicit none
   type(mpas_pool_type), pointer, intent(inout) :: pool

   if (associated(pool)) then
      call mpas_pool_destroy_pool(pool)
   end if

end subroutine delete_pool

! ------------------------------------------------------------------------------

subroutine copy_fields(self,rhs)

   implicit none
   class(mpas_fields), intent(inout) :: self
   class(mpas_fields), intent(in)    :: rhs
   type (MPAS_Time_type) :: rhs_time
   integer :: ierr

!   write(*,*)'--> copy_fields: copy subFields Pool'

   self % nf = rhs % nf
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   allocate(self % fldnames(self % nf))
   self % fldnames(:) = rhs % fldnames(:)

   self % nf_ci = rhs % nf_ci
   if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)
   allocate(self % fldnames_ci(self % nf_ci))
   self % fldnames_ci(:) = rhs % fldnames_ci(:)

   rhs_time = mpas_get_clock_time(rhs % clock, MPAS_NOW, ierr)
   call mpas_set_clock_time(self % clock, rhs_time, MPAS_NOW)

   call copy_pool(rhs % subFields, self % subFields)

!   write(*,*)'--> copy_fields done'

end subroutine copy_fields

! ------------------------------------------------------------------------------

subroutine copy_pool(pool_src, pool)

   implicit none
   type(mpas_pool_type), pointer, intent(in)    :: pool_src
   type(mpas_pool_type), pointer, intent(inout) :: pool

   ! Duplicate the members of pool_src into pool and
   ! do a deep copy of the fields
   call delete_pool(pool)
   call mpas_pool_create_pool(pool)
   call mpas_pool_clone_pool(pool_src, pool)

end subroutine copy_pool

! ------------------------------------------------------------------------------

subroutine read_fields(self, f_conf, vdate)

   implicit none
   class(mpas_fields),        intent(inout) :: self     !< Field
   type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
   type(datetime),            intent(inout) :: vdate    !< DateTime
   character(len=:), allocatable :: str
   character(len=20)       :: sdate
   type (MPAS_Time_type)   :: local_time
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   integer                 :: ierr = 0, ngrid
   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer    :: pressure, pressure_base, pressure_p

!   write(*,*)'--> read_fields'
   call f_conf%get_or_die("date",str)
   sdate = str
   call datetime_set(sdate, vdate)

   call f_conf%get_or_die("filename",str)
   call swap_name_member(f_conf, str)
   temp_filename = str
!   write(*,*)'--> read_fields: Reading ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   ! TODO: we can get streamID from yaml
   !streamID = 'restart'
   !streamID = 'input'
   streamID = 'output'
   !streamID = 'da'
   ierr = 0
   self % manager => self % geom % domain % streamManager
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(sdate,dateTimeString,1)
!   write(*,*)'--> read_fields: dateTimeString: ',trim(dateTimeString)
   call mpas_set_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
   call mpas_set_clock_time(self % clock, local_time, MPAS_NOW)
   call mpas_set_clock_time(self % geom % domain % clock, local_time, MPAS_START_TIME)
   call mpas_expand_string(dateTimeString, -1, temp_filename, filename)
   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
!   write(*,*)'--> read_fields: Reading ',trim(filename)
   call MPAS_stream_mgr_read(self % manager, streamID=streamID, &
                           & when=dateTimeString, rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      write(message,*) '--> read_fields: MPAS_stream_mgr_read failed ierr=',ierr
      call abor1_ftn(message)
   end if

   !==TODO: Speific part when reading parameterEst. for BUMP.	
   !      : They write/read a list of variables directly.	
   if (f_conf%has("no_transf")) then
      call f_conf%get_or_die("no_transf",ierr)
      if(ierr .eq. 1) then
         call da_copy_all2sub_fields(self % geom % domain, self % subFields)
        return
      endif
   endif

   !(1) diagnose pressure
   call mpas_pool_get_subpool(self % geom % domain % blocklist % structs, 'diag', diag)
   call mpas_pool_get_field(diag, 'pressure_p', pressure_p)
   call mpas_pool_get_field(diag, 'pressure_base', pressure_base)
   call mpas_pool_get_field(diag, 'pressure', pressure)
   ngrid = self % geom % nCellsSolve
   pressure%array(:,1:ngrid) = pressure_base%array(:,1:ngrid) + pressure_p%array(:,1:ngrid)

   !(2) copy all to subFields & diagnose temperature
   call update_diagnostic_fields(self % geom % domain, self % subFields, self % geom % nCellsSolve)

end subroutine read_fields


subroutine update_diagnostic_fields(domain, subFields, ngrid)

   implicit none
   type (domain_type), pointer,    intent(inout) :: domain
   type (mpas_pool_type), pointer, intent(inout) :: subFields
   integer,                        intent(in)    :: ngrid
   type (field2DReal), pointer    :: theta, pressure, temperature, specific_humidity
   type (field3DReal), pointer    :: scalars
   type (mpas_pool_type), pointer :: state
   integer, pointer :: index_qv

   !(1) copy all to subFields
   call da_copy_all2sub_fields(domain, subFields)

   !(2) diagnose temperature
   !Special case: Convert theta and pressure to temperature
   !              Convert water vapor mixing ratio to specific humidity [ q = w / (1 + w) ]
   !NOTE: This formula is somewhat different with MPAS one's (in physics, they use "exner")
   !    : If T diagnostic is added in, for example, subroutine atm_compute_output_diagnostics,
   !    : we need to include "exner" in stream_list.for.reading

   call mpas_pool_get_field(domain % blocklist % allFields, 'theta', theta)
   call mpas_pool_get_field(domain % blocklist % allFields, 'pressure', pressure)
   call mpas_pool_get_field(subFields, 'temperature', temperature)
   call mpas_pool_get_field(domain % blocklist % allFields, 'scalars', scalars)
   call mpas_pool_get_field(subFields, 'spechum', specific_humidity)

   call mpas_pool_get_subpool(domain % blocklist % structs,'state',state)
   call mpas_pool_get_dimension(state, 'index_qv', index_qv)

   call theta_to_temp(theta % array(:,1:ngrid), pressure % array(:,1:ngrid), temperature % array(:,1:ngrid))
   call w_to_q( scalars % array(index_qv,:,1:ngrid) , specific_humidity % array(:,1:ngrid) )

end subroutine update_diagnostic_fields

! ------------------------------------------------------------------------------

subroutine write_fields(self, f_conf, vdate)

   implicit none
   class(mpas_fields),        intent(inout) :: self   !< Field
   type(fckit_configuration), intent(in)    :: f_conf !< Configuration
   type(datetime),            intent(in)    :: vdate  !< DateTime
   character(len=:), allocatable :: str
   character(len=20)       :: validitydate
   integer                 :: ierr
   type (MPAS_Time_type)   :: fld_time, write_time
   character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename, temp_filename

   call da_copy_sub2all_fields(self % geom % domain, self % subFields)

   call datetime_to_string(vdate, validitydate)
!   write(*,*)'--> write_fields: ',trim(validitydate)
   call f_conf%get_or_die("filename",str)
   call swap_name_member(f_conf, str)
   temp_filename = str
!   write(*,*)'--> write_fields: ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(validitydate,dateTimeString,-1)
   ierr = 0
   call mpas_set_time(write_time, dateTimeString=dateTimeString, ierr=ierr)
   fld_time = mpas_get_clock_time(self % clock, MPAS_NOW, ierr)
   call mpas_get_time(fld_time, dateTimeString=dateTimeString2, ierr=ierr)
   ! write(*,*)'check time --> write_fields: write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)

   self % manager => self % geom % domain % streamManager
   ! TODO: we can get streamID from yaml
   ! TODO: should we pick different stream lists for mpas_state and mpas_increment?
   !streamID = 'restart'
   streamID = 'output'
   !streamID = 'da'
   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)

!   write(*,*)'--> write_fields: writing ',trim(filename)
   call mpas_stream_mgr_write(self % geom % domain % streamManager, streamID=streamID, &
        forceWriteNow=.true., writeTime=dateTimeString, ierr=ierr)
   if ( ierr .ne. 0  ) then
     write(message,*) '--> write_fields: MPAS_stream_mgr_write failed ierr=',ierr
     call abor1_ftn(message)
   end if

end subroutine write_fields

! ------------------------------------------------------------------------------

subroutine change_resol_fields(self,rhs)

   implicit none
   class(mpas_fields), intent(inout) :: self
   class(mpas_fields), intent(in)    :: rhs

   if (self%geom%nCells == rhs%geom%nCells .and.  self%geom%nVertLevels == rhs%geom%nVertLevels) then
     call self%copy(rhs)
   else if (self%geom%nVertLevels == rhs%geom%nVertLevels) then
      call interpolate_fields(self,rhs)
   else
     write(message,*) '--> change_resol_fields: ',self%geom%nCells, rhs%geom%nCells, self%geom%nVertLevels, rhs%geom%nVertLevels
     call fckit_log%info(message)
     call abor1_ftn("mpas_fields_mod:change_resol_fields: VertLevels dimension mismatch")
   endif

end subroutine change_resol_fields

! ------------------------------------------------------------------------------

subroutine zeros_(self)

   implicit none
   class(mpas_fields), intent(inout) :: self

   call da_constant(self % subFields, MPAS_JEDI_ZERO_kr, fld_select = self % fldnames_ci)

end subroutine zeros_

! ------------------------------------------------------------------------------

subroutine ones_(self)

   implicit none
   class(mpas_fields), intent(inout) :: self

   call da_constant(self % subFields, MPAS_JEDI_ONE_kr, fld_select = self % fldnames_ci)

end subroutine ones_

! ------------------------------------------------------------------------------

subroutine random_(self)

   implicit none
   class(mpas_fields), intent(inout) :: self

   call da_random(self % subFields, fld_select = self % fldnames_ci)

end subroutine random_

! ------------------------------------------------------------------------------

subroutine gpnorm_(self, nf, pstat)

   implicit none
   class(mpas_fields),   intent(in)  :: self
   integer,              intent(in)  :: nf
   real(kind=kind_real), intent(out) :: pstat(3, nf)

   call da_gpnorm(self % subFields, self % geom % domain % dminfo, nf, pstat, fld_select = self % fldnames_ci(1:nf))

end subroutine gpnorm_

! ------------------------------------------------------------------------------

subroutine rms_(self, prms)

   implicit none
   class(mpas_fields),   intent(in)  :: self
   real(kind=kind_real), intent(out) :: prms

   call da_fldrms(self % subFields, self % geom % domain % dminfo, prms, fld_select = self % fldnames_ci)

end subroutine rms_

! ------------------------------------------------------------------------------

subroutine self_add_(self,rhs)

   implicit none
   class(mpas_fields), intent(inout) :: self
   class(mpas_fields), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'add'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_add_

! ------------------------------------------------------------------------------

subroutine self_schur_(self,rhs)

   implicit none
   class(mpas_fields), intent(inout) :: self
   class(mpas_fields), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'schur'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_schur_

! ------------------------------------------------------------------------------

subroutine self_sub_(self,rhs)

   implicit none
   class(mpas_fields), intent(inout) :: self
   class(mpas_fields), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'sub'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_sub_

! ------------------------------------------------------------------------------

subroutine self_mult_(self,zz)

   implicit none
   class(mpas_fields),   intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz

   call da_self_mult(self % subFields, zz)

end subroutine self_mult_

! ------------------------------------------------------------------------------

subroutine axpy_(self,zz,rhs)

   implicit none
   class(mpas_fields),   intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz
   class(mpas_fields),   intent(in)    :: rhs

   call da_axpy(self % subFields, rhs % subFields, zz, fld_select = self % fldnames_ci)

end subroutine axpy_

! ------------------------------------------------------------------------------

subroutine dot_prod_(self,fld,zprod)

   implicit none
   class(mpas_fields),    intent(in)    :: self, fld
   real(kind=kind_real),  intent(inout) :: zprod

   call da_dot_product(self % subFields, fld % subFields, self % geom % domain % dminfo, zprod)

end subroutine dot_prod_

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)

   implicit none
   character(len=2),     intent(in) :: cop
   class(mpas_fields),   intent(in) :: fld
   type(ufo_locations),  intent(in) :: locs
   type(oops_variables), intent(in) :: vars
   type(ufo_geovals),    intent(in) :: gom
   integer :: jvar
   character(len=26) :: cinfo

   cinfo="mpas_fields:checks "//cop//" : "

   !Check things are the sizes we expect
   !------------------------------------
   if( gom%nvar .ne. vars%nvars() )then
      call abor1_ftn(cinfo//"nvar wrong size")
   endif
   if( .not. allocated(gom%geovals) )then
      call abor1_ftn(cinfo//"geovals unallocated")
   endif
   if( size(gom%geovals) .ne. vars%nvars() )then
      call abor1_ftn(cinfo//"geovals wrong size")
   endif
   if (cop/="tl" .and. cop/='nl') then !BJJ or cop='ad'.   why only check ad??? because ad is set/allocated in UFO side ??
                                                           ! also, the dimension is not general for all possible geovals.
   if (.not.gom%linit) then
      call abor1_ftn(cinfo//"geovals not initialized")
   endif
   if (.not.allocated(gom%geovals(jvar)%vals)) then
     call abor1_ftn(cinfo//"vals not allocated")
   endif

   endif

!   write(*,*)'--> interp_checks: ',cinfo,' done'

end subroutine interp_checks

! ------------------------------------------------------------------------------

!> \brief Populates subfields of self using rhs
!!
!! \details **interpolate_fields** This subroutine is called when creating
!! a new mpas_fields type (self) using an existing mpas_fields (rhs) as a source that
!! has a different geometry/mesh (but the same number of VertLevels). It populates
!! the subfields of self, interpolating the data from rhs. It can use either
!! bump or unstructured interpolation for the interpolation routine.
subroutine interpolate_fields(self,rhs)

  implicit none
  class(mpas_fields), intent(inout) :: self !< mpas_fields being populated
  class(mpas_fields), intent(in)    :: rhs  !< mpas_fields used as source

  type(bump_interpolator) :: bumpinterp
  type(unstrc_interp)     :: unsinterp
  type (mpas_pool_iterator_type) :: poolItr
  real(kind=kind_real), allocatable :: interp_in(:,:), interp_out(:,:)
  real (kind=kind_real), dimension(:), pointer :: r1d_ptr
  real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr
  integer, dimension(:), pointer :: i1d_ptr
  integer, dimension(:,:), pointer :: i2d_ptr
  integer :: rhs_nCells, self_nCells, maxlevels, nlevels, jlev
  logical             :: use_bump_interp
  integer, allocatable :: rhsDims(:)

  use_bump_interp = rhs%geom%use_bump_interpolation

  if (use_bump_interp) then
    call initialize_bumpinterp(self%geom, rhs%geom, bumpinterp)
  else
    call initialize_uns_interp(self%geom, rhs%geom, unsinterp)
  endif

  self % nf = rhs % nf
  if (allocated(self % fldnames)) deallocate(self % fldnames)
  allocate(self % fldnames(self % nf))
  self % fldnames(:) = rhs % fldnames(:)

  self % nf_ci = rhs % nf_ci
  if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)
  allocate(self % fldnames_ci(self % nf_ci))
  self % fldnames_ci(:) = rhs % fldnames_ci(:)

  call create_pool(self % geom % domain, self % nf, self % fldnames, self % subFields)

  ! Interpolate field from rhs mesh to self mesh using pre-calculated weights
  ! ----------------------------------------------------------------------------------
  maxlevels = rhs%geom%nVertLevelsP1
  rhs_nCells = rhs%geom%nCellsSolve
  self_nCells = self%geom%nCellsSolve

  allocate(interp_in(rhs_nCells, maxlevels))
  allocate(interp_out(self_nCells, maxlevels))

  call mpas_pool_begin_iteration(rhs%subFields)
  do while ( mpas_pool_get_next_member(rhs%subFields, poolItr) )
  if (poolItr % memberType == MPAS_POOL_FIELD) then
    write(message,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
    call fckit_log%debug(message)
    ! Get the rhs fields out of the rhs%subFields pool and into an array that
    ! can be passed to the interpolator apply functions.
    ! (Four cases to cover: 1D/2D and real/integer variable.)
    if (poolItr % nDims == 1) then
      nlevels = 1
      if (poolItr % dataType == MPAS_POOL_INTEGER) then
        call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), i1d_ptr)
        interp_in(:,1) = real( i1d_ptr(1:rhs_nCells), kind_real)
      else if (poolItr % dataType == MPAS_POOL_REAL) then
        call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), r1d_ptr)
        interp_in(:,1) = r1d_ptr(1:rhs_nCells)
      endif
    else if (poolItr % nDims == 2) then
      rhsDims = getSolveDimensions(rhs%subFields, poolItr)
      nlevels = rhsDims(1)
      if (nlevels > maxlevels) then
        write(message,*) '--> interpolate_fields: nlevels > maxlevels, ', nlevels, maxlevels
        call abor1_ftn(message)
      endif
      if (poolItr % dataType == MPAS_POOL_INTEGER) then
        call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), i2d_ptr)
        interp_in(1:rhs_nCells,1:nlevels) = real( transpose (i2d_ptr(1:nlevels,1:rhs_nCells)), kind_real )
      else if (poolItr % dataType == MPAS_POOL_REAL) then
        call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), r2d_ptr)
        interp_in(1:rhs_nCells,1:nlevels) = transpose(r2d_ptr(1:nlevels,1:rhs_nCells))
      endif
    else
      write(message,*) '--> interpolate_fields: poolItr % nDims == ',poolItr % nDims,' not handled'
      call abor1_ftn(message)
    endif

    if (use_bump_interp) then
      call bumpinterp%apply(interp_in(:,1:nlevels), &
                            interp_out(:,1:nlevels), &
                            trans_in=.false.)
    else
      do jlev = 1, nlevels
        call unsinterp%apply(interp_in(:,jlev),interp_out(:,jlev))
      end do
    endif

    ! Put the interpolated results into the self%subFields pool
    if (poolItr % nDims == 1) then
      if (poolItr % dataType == MPAS_POOL_INTEGER) then
        call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i1d_ptr)
        i1d_ptr(1:self_nCells) = int (interp_out(:,1))
      else if (poolItr % dataType == MPAS_POOL_REAL) then
        call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r1d_ptr)
        r1d_ptr(1:self_nCells) = interp_out(:,1)
      endif
    else if (poolItr % nDims == 2) then
      if (poolItr % dataType == MPAS_POOL_INTEGER) then
        call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i2d_ptr)
        i2d_ptr(1:nlevels,1:self_nCells) = transpose (int (interp_out(1:self_nCells,1:nlevels)))
      else if (poolItr % dataType == MPAS_POOL_REAL) then
        call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r2d_ptr)
        r2d_ptr(1:nlevels,1:self_nCells) = transpose (interp_out(1:self_nCells,1:nlevels))
      endif
    endif
  endif
  end do !- end of pool iteration
  deallocate(interp_in)
  deallocate(interp_out)
  if (allocated(rhsDims)) deallocate(rhsDims)

end subroutine interpolate_fields
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> \brief Initializes a bump interpolation type
!!
!! \details **initialize_bumpinterp** This subroutine calls bumpinterp%init,
!! which calculates the barycentric weights used to interpolate data between the
!! geom_from locations and the geom_to locations.
subroutine initialize_bumpinterp(geom_to, geom_from, bumpinterp)

   implicit none
   class(mpas_geom), intent(in)           :: geom_to     !< geometry interpolating to
   class(mpas_geom), intent(in)           :: geom_from   !< geometry interpolating from
   type(bump_interpolator), intent(inout) :: bumpinterp  !< bump interpolator

   real(kind=kind_real), allocatable :: lats_to(:), lons_to(:)

   allocate( lats_to(geom_to%nCellsSolve) )
   allocate( lons_to(geom_to%nCellsSolve) )
   lats_to(:) = geom_to%latCell( 1:geom_to%nCellsSolve ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
   lons_to(:) = geom_to%lonCell( 1:geom_to%nCellsSolve ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees

   call bumpinterp%init(geom_from%f_comm,afunctionspace_in=geom_from%afunctionspace,lon_out=lons_to,lat_out=lats_to, &
      & nl=geom_from%nVertLevels)

   ! Release memory
   ! --------------
   deallocate(lats_to)
   deallocate(lons_to)

end subroutine initialize_bumpinterp

! --------------------------------------------------------------------------------------------------

!> \brief Initializes an unstructured interpolation type
!!
!! \details **initialize_uns_interp** This subroutine calls unsinterp%create,
!! which calculates the barycentric weights used to interpolate data between the
!! geom_from locations and the geom_to locations.
subroutine initialize_uns_interp(geom_to, geom_from, unsinterp)

   implicit none
   class(mpas_geom), intent(in)           :: geom_to     !< geometry interpolating to
   class(mpas_geom), intent(in)           :: geom_from   !< geometry interpolating from
   type(unstrc_interp),     intent(inout) :: unsinterp   !< unstructured interpolator

   integer :: nn, ngrid_from, ngrid_to
   character(len=8) :: wtype = 'barycent'
   real(kind=kind_real), allocatable :: lats_from(:), lons_from(:), lats_to(:), lons_to(:)

   ! Get the Solution dimensions
   ! ---------------------------
   ngrid_from = geom_from%nCellsSolve
   ngrid_to   = geom_to%nCellsSolve

   !Calculate interpolation weight
   !------------------------------------------
   allocate( lats_from(ngrid_from) )
   allocate( lons_from(ngrid_from) )
   lats_from(:) = geom_from%latCell( 1:ngrid_from ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
   lons_from(:) = geom_from%lonCell( 1:ngrid_from ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
   allocate( lats_to(ngrid_to) )
   allocate( lons_to(ngrid_to) )
   lats_to(:) = geom_to%latCell( 1:ngrid_to ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
   lons_to(:) = geom_to%lonCell( 1:ngrid_to ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees

   ! Initialize unsinterp
   ! ---------------
   nn = 3 ! number of nearest neigbors
   call unsinterp%create(geom_from%f_comm, nn, wtype, &
                         ngrid_from, lats_from, lons_from, &
                         ngrid_to, lats_to, lons_to)

   ! Release memory
   ! --------------
   deallocate(lats_from)
   deallocate(lons_from)
   deallocate(lats_to)
   deallocate(lons_to)

end subroutine initialize_uns_interp

subroutine serial_size(self, vsize)

   implicit none

   ! Passed variables
   class(mpas_fields),intent(in) :: self
   integer,intent(out) :: vsize !< Size

   ! Local variables
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: solveDims(:)

   ! Initialize
   vsize = 0

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         solveDims = getSolveDimensions(self%subFields, poolItr)
         vsize = vsize + product(solveDims)
         deallocate(solveDims)
      endif
   enddo

end subroutine serial_size

! ------------------------------------------------------------------------------

subroutine serialize_fields(self, vsize, vect_inc)

   implicit none

   ! Passed variables
   class(mpas_fields),intent(in) :: self             !< Increment
   integer,intent(in) :: vsize                      !< Size
   real(kind_real),intent(out) :: vect_inc(vsize)   !< Vector

   ! Local variables
   integer :: index, nvert, nhoriz, vv, hh
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: solveDims(:)

   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
   integer, dimension(:), pointer :: i1d_ptr_a
   integer, dimension(:,:), pointer :: i2d_ptr_a

   ! Initialize
   index = 0

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         solveDims = getSolveDimensions(self%subFields, poolItr)
         if (poolItr % nDims == 1) then
            nhoriz = solveDims(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i1d_ptr_a)
               do hh = 1, nhoriz
                  vect_inc(index + 1) = real(i1d_ptr_a(hh), kind=kind_real)
                  index = index + 1
               enddo
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r1d_ptr_a)
               do hh = 1, nhoriz
                  vect_inc(index + 1) = r1d_ptr_a(hh)
                  index = index + 1
               enddo
            endif
         elseif (poolItr % nDims == 2) then
            nvert = solveDims(1)
            nhoriz = solveDims(2)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     vect_inc(index + 1) = real(i2d_ptr_a(vv, hh), kind=kind_real)
                     index = index + 1
                  enddo
               enddo
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     vect_inc(index + 1) = r2d_ptr_a(vv, hh)
                     index = index + 1
                  enddo
               enddo
            endif
         else
            write(message,*) '--> serialize_fields: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         deallocate(solveDims)
      endif
   enddo

end subroutine serialize_fields

! --------------------------------------------------------------------------------------------------

subroutine deserialize_fields(self, vsize, vect_inc, index)

   implicit none

   ! Passed variables
   class(mpas_fields),intent(inout) :: self               !< Increment
   integer,intent(in) :: vsize                           !< Size
   real(kind_real),intent(in) :: vect_inc(vsize)         !< Vector
   integer,intent(inout) :: index                        !< Index

   ! Local variables
   integer :: nvert, nhoriz, vv, hh
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: solveDims(:)

   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
   integer, dimension(:), pointer :: i1d_ptr_a
   integer, dimension(:,:), pointer :: i2d_ptr_a

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         solveDims = getSolveDimensions(self%subFields, poolItr)
         if (poolItr % nDims == 1) then
            nhoriz = solveDims(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i1d_ptr_a)
               do hh = 1, nhoriz
                  i1d_ptr_a(hh) = int ( vect_inc(index + 1) )
                  index = index + 1
               enddo
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r1d_ptr_a)
               do hh = 1, nhoriz
                  r1d_ptr_a(hh) = vect_inc(index + 1)
                  index = index + 1
               enddo
            endif
         elseif (poolItr % nDims == 2) then
            nvert = solveDims(1)
            nhoriz = solveDims(2)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     i2d_ptr_a(vv, hh) = int ( vect_inc(index + 1) )
                     index = index + 1
                  enddo
               enddo
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     r2d_ptr_a(vv, hh) = vect_inc(index + 1)
                     index = index + 1
                  enddo
               enddo
            endif
         else
            write(message,*) '--> deserialize_fields: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         deallocate(solveDims)
      endif
   enddo

end subroutine deserialize_fields

function has_field(self, fieldname) result(has)
   class(mpas_fields), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   logical :: has
   has = (ufo_vars_getindex(self % fldnames, fieldname) > 0)
end function has_field

function has_fields(self, fieldnames) result(has)
   class(mpas_fields), intent(in) :: self
   character(len=*), intent(in) :: fieldnames(:)
   integer :: i
   logical :: has
   has = .true.
   do i = 1, size(fieldnames)
      has = has .and. self%has(fieldnames(i))
   end do
end function has_fields

end module mpas_fields_mod
