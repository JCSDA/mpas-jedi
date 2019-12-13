! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_field_utils_mod

use fckit_configuration_module, only: fckit_configuration

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables

!ufo
use ufo_vars_mod, only: MAXVARLEN

!MPAS-Model
use atm_core, only: atm_simulation_clock_init, atm_compute_output_diagnostics
use mpas_constants
use mpas_derived_types
use mpas_kind_types, only: StrKIND
use mpas_pool_routines
use mpas_stream_manager
use mpas_timekeeping

!mpas-jedi
use mpas_geom_mod
use mpas4da_mod


private

public :: mpas_field, mpas_field_registry, &
          interp_checks, &
          create_field, delete_field, &
          copy_field, copy_pool, &
          update_diagnostic_fields, &
          mpas_hydrometeor_fields,  &
          mpas_re_fields

! ------------------------------------------------------------------------------

   !> Fortran derived type to hold MPAS field
   type :: mpas_field
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

     procedure :: change_resol => change_resol_field
     procedure :: copy         => copy_field
     procedure :: create       => create_field
     procedure :: delete       => delete_field
     procedure :: read_file    => read_field
     procedure :: write_file   => write_field

   end type mpas_field

!   abstract interface
!
!   ! ------------------------------------------------------------------------------
!
!      subroutine read_file_(self, f_conf, vdate)
!         import mpas_field, fckit_configuration, datetime
!         implicit none
!         class(mpas_field),         intent(inout) :: self
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

subroutine create_field(self, geom, vars, vars_ci)

    implicit none

    class(mpas_field), intent(inout)       :: self
    type(mpas_geom),   intent(in), pointer :: geom
    type(oops_variables), intent(in)       :: vars, vars_ci

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

!    write(*,*)'--> create_field: self % nf =',self % nf
!    write(*,*)'--> create_field: allocate ::',self % fldnames(:)
 
    ! link geom
    if (associated(geom)) then
      self % geom => geom
    else
      write(*,*)'--> create_field: geom not associated'
      call abor1_ftn("--> create_field: geom not associated")
    end if

    ! clock creation
    allocate(self % clock)
    call atm_simulation_clock_init(self % clock, self % geom % domain % blocklist % configs, ierr)
    if ( ierr .ne. 0 ) then
       call abor1_ftn("--> create_field: atm_simulation_clock_init problem")
    end if

!    write(*,*)'--> create_field: sub Pool from list of variable ',self % nf
    call create_pool(self % geom % domain, self % nf, self % fldnames, self % subFields)

    call self%zeros() !-- set zero for self % subFields

    return

end subroutine create_field

! ------------------------------------------------------------------------------

subroutine create_pool(domain, nf_in, fldnames, pool)

    implicit none
    type (domain_type), pointer, intent(in) :: domain
    integer, intent(in) :: nf_in
    character (len=*), intent(in) :: fldnames(:)
    type (mpas_pool_type), pointer, intent(out) :: pool

    integer :: nfields
    character(len=1024) :: buf

    call da_make_subpool(domain, pool, nf_in, fldnames, nfields)

    if ( nf_in .ne. nfields  ) then
       write(buf,*) "--> create_pool: dimension mismatch ", nf_in, nfields
       call abor1_ftn(buf)
    end  if

end subroutine create_pool

! ------------------------------------------------------------------------------

subroutine delete_field(self)

   implicit none
   class(mpas_field), intent(inout) :: self
   integer :: ierr = 0 
 
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)

!   write(*,*)'--> delete_field: deallocate subFields Pool'
   call delete_pool(self % subFields)

   call mpas_destroy_clock(self % clock, ierr)
   if ( ierr .ne. 0  ) then
      write(*,*) '--> delete_field deallocate clock failed'
   end if
!   write(*,*)'--> delete_field done'

   return

end subroutine delete_field

! ------------------------------------------------------------------------------

subroutine delete_pool(pool)

   implicit none
   type(mpas_pool_type), pointer, intent(inout) :: pool
  
   if (associated(pool)) then
      call mpas_pool_destroy_pool(pool)
   end if

end subroutine delete_pool

! ------------------------------------------------------------------------------

subroutine copy_field(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self
   class(mpas_field), intent(in)    :: rhs

!   write(*,*)'--> copy_field: copy subFields Pool'

   self % nf = rhs % nf
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   allocate(self % fldnames(self % nf))
   self % fldnames(:) = rhs % fldnames(:)

   self % nf_ci = rhs % nf_ci
   if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)
   allocate(self % fldnames_ci(self % nf_ci))
   self % fldnames_ci(:) = rhs % fldnames_ci(:)

   call copy_pool(rhs % subFields, self % subFields)

!   write(*,*)'--> copy_field done'
  
end subroutine copy_field

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

subroutine read_field(self, f_conf, vdate)

   implicit none
   class(mpas_field),         intent(inout) :: self     !< Field
   type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
   type(datetime),            intent(inout) :: vdate    !< DateTime
   character(len=:), allocatable :: str
   character(len=20)       :: sdate
   type (MPAS_Time_type)   :: local_time
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   integer                 :: ierr = 0, ngrid
   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer    :: pressure, pressure_base, pressure_p
   character(len=1024) :: buf

!   write(*,*)'--> read_field'
   call f_conf%get_or_die("date",str)
   sdate = str
   call datetime_set(sdate, vdate)

   call f_conf%get_or_die("filename",str)
   temp_filename = str
!   write(*,*)'--> read_field: Reading ',trim(temp_filename)
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
!   write(*,*)'--> read_field: dateTimeString: ',trim(dateTimeString)
   call mpas_set_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
   !call mpas_set_clock_time(self % clock, local_time, MPAS_START_TIME)
   call mpas_set_clock_time(self % geom % domain % clock, local_time, MPAS_START_TIME)
   call mpas_expand_string(dateTimeString, -1, temp_filename, filename)
   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(*,*)'--> read_field: Reading ',trim(filename)
   call MPAS_stream_mgr_read(self % manager, streamID=streamID, &
                           & when=dateTimeString, rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      call abor1_ftn(buf)
      write(buf,*) '--> read_field: MPAS_stream_mgr_read failed ierr=',ierr
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

end subroutine read_field


subroutine update_diagnostic_fields(domain, subFields, ngrid)

   use mpas2ufo_vars_mod, only: w_to_q, theta_to_temp

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
   call mpas_pool_get_field(subFields, 'pressure', pressure)
   call mpas_pool_get_field(subFields, 'temperature', temperature)
   call mpas_pool_get_field(domain % blocklist % allFields, 'scalars', scalars)
   call mpas_pool_get_field(subFields, 'spechum', specific_humidity)

   call mpas_pool_get_subpool(domain % blocklist % structs,'state',state)
   call mpas_pool_get_dimension(state, 'index_qv', index_qv)

   call theta_to_temp(theta % array(:,1:ngrid), pressure % array(:,1:ngrid), temperature % array(:,1:ngrid))
   call w_to_q( scalars % array(index_qv,:,1:ngrid) , specific_humidity % array(:,1:ngrid) )

end subroutine update_diagnostic_fields

! ------------------------------------------------------------------------------

subroutine write_field(self, f_conf, vdate)

   implicit none
   class(mpas_field),         intent(inout) :: self   !< Field
   type(fckit_configuration), intent(in)    :: f_conf !< Configuration
   type(datetime),            intent(in)    :: vdate  !< DateTime
   character(len=:), allocatable :: str
   character(len=20)       :: validitydate
   integer                 :: ierr, iskip
   type (MPAS_Time_type)   :: fld_time, write_time
   character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename, temp_filename
   character(len=1024)     :: buf

   call da_copy_sub2all_fields(self % geom % domain, self % subFields)

   call datetime_to_string(vdate, validitydate)
!   write(*,*)'--> write_field: ',trim(validitydate)
   call f_conf%get_or_die("filename",str)
   temp_filename = str
!   write(*,*)'--> write_field: ',trim(temp_filename)
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
   iskip = 0
   if (f_conf%has("SkipMPASTimeCheck")) then
     call f_conf%get_or_die("SkipMPASTimeCheck",iskip)
   end if
   if(iskip.ne.1) then
     if ( fld_time .NE. write_time ) then
        write(*,*)'--> write_field: write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
        call abor1_ftn('Different times MPAS_stream_mgr_write failed ')
     end if
   end if
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)
   ! Filename for writing out ensembles
   if (f_conf%has("member")) then
     call f_conf%get_or_die("member",str)
     temp_filename = str
     filename=trim(filename) // '.' // trim(temp_filename)
   endif
   self % manager => self % geom % domain % streamManager
   ! TODO: we can get streamID from yaml
   ! TODO: should we pick different stream lists for mpas_state and mpas_increment?
   !streamID = 'restart'
   streamID = 'output' 
   !streamID = 'da'
   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)

   write(*,*)'--> write_field: writing ',trim(filename)
   call mpas_stream_mgr_write(self % geom % domain % streamManager, streamID=streamID, forceWriteNow=.true., ierr=ierr)
   if ( ierr .ne. 0  ) then
     write(buf,*) '--> write_field: MPAS_stream_mgr_write failed ierr=',ierr
     call abor1_ftn(buf)
   end if

end subroutine write_field

! ------------------------------------------------------------------------------

subroutine change_resol_field(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self
   class(mpas_field), intent(in)    :: rhs

   ! FIXME: We just copy rhs to self for now. Need an actual interpolation routine later. (SH)
   if (self%geom%nCells == rhs%geom%nCells .and.  self%geom%nVertLevels == rhs%geom%nVertLevels) then
     call self%copy(rhs)
   else
     write(*,*) '--> write_field: ',self%geom%nCells, rhs%geom%nCells, self%geom%nVertLevels, rhs%geom%nVertLevels
     call abor1_ftn("change_resol_field: dimension mismatch")
   endif

end subroutine change_resol_field

! ------------------------------------------------------------------------------

subroutine zeros_(self)

   implicit none
   class(mpas_field), intent(inout) :: self

   call da_zeros(self % subFields, fld_select = self % fldnames_ci)

end subroutine zeros_

! ------------------------------------------------------------------------------

subroutine random_(self)

   implicit none
   class(mpas_field), intent(inout) :: self
   
   call da_random(self % subFields, fld_select = self % fldnames_ci)

end subroutine random_

! ------------------------------------------------------------------------------

subroutine gpnorm_(self, nf, pstat)

   implicit none
   class(mpas_field),    intent(in)  :: self
   integer,              intent(in)  :: nf
   real(kind=kind_real), intent(out) :: pstat(3, nf)

   call da_gpnorm(self % subFields, self % geom % domain % dminfo, nf, pstat, fld_select = self % fldnames_ci(1:nf))

end subroutine gpnorm_

! ------------------------------------------------------------------------------

subroutine rms_(self, prms)

   implicit none
   class(mpas_field),    intent(in)  :: self
   real(kind=kind_real), intent(out) :: prms

   call da_fldrms(self % subFields, self % geom % domain % dminfo, prms, fld_select = self % fldnames_ci)

end subroutine rms_

! ------------------------------------------------------------------------------

subroutine self_add_(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self
   class(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'add'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_add_

! ------------------------------------------------------------------------------

subroutine self_schur_(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self
   class(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'schur'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_schur_

! ------------------------------------------------------------------------------

subroutine self_sub_(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self
   class(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'sub'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = self % fldnames_ci)

end subroutine self_sub_

! ------------------------------------------------------------------------------

subroutine self_mult_(self,zz)

   implicit none
   class(mpas_field),    intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz

   call da_self_mult(self % subFields, zz)

end subroutine self_mult_

! ------------------------------------------------------------------------------

subroutine axpy_(self,zz,rhs)

   implicit none
   class(mpas_field),    intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz
   class(mpas_field),    intent(in)    :: rhs

   call da_axpy(self % subFields, rhs % subFields, zz, fld_select = self % fldnames_ci)

end subroutine axpy_

! ------------------------------------------------------------------------------

subroutine dot_prod_(self,fld,zprod)

   implicit none
   class(mpas_field),     intent(in)    :: self, fld
   real(kind=kind_real),  intent(inout) :: zprod

   call da_dot_product(self % subFields, fld % subFields, self % geom % domain % dminfo, zprod)

end subroutine dot_prod_

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)

   use ufo_locs_mod, only: ufo_locs
   use ufo_geovals_mod, only: ufo_geovals
   implicit none
   character(len=2),     intent(in) :: cop
   class(mpas_field),    intent(in) :: fld
   type(ufo_locs),       intent(in) :: locs
   type(oops_variables), intent(in) :: vars
   type(ufo_geovals),    intent(in) :: gom
   integer :: jvar
   character(len=26) :: cinfo

   cinfo="mpas_field:checks "//cop//" : "

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

end module mpas_field_utils_mod
