! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpas_fields_mod

! atlas
use atlas_module, only: atlas_field, atlas_fieldset, atlas_real, atlas_metadata

! fckit
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use string_utils, only: swap_name_member
use slow_oops_unstructured_interpolation_mod, only : unstrc_interp

! saber
use interpolatorbump_mod, only: bump_interpolator

!ufo
use ufo_vars_mod, only: MAXVARLEN, ufo_vars_getindex, var_prsi
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
use mpas_dmpar, only : mpas_dmpar_exch_halo_field, mpas_dmpar_exch_halo_adj_field

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas4da_mod
use mpas2ufo_vars_mod, only: w_to_q, theta_to_temp
use mpas_kinds, only : c_real_type

implicit none

private

public :: mpas_fields, mpas_fields_registry, &
          create_fields, delete_fields, &
          copy_fields, copy_pool, &
          update_diagnostic_fields, &
          mpas_hydrometeor_fields,  &
          mpas_re_fields, &
          cellCenteredWindFields, &
          moistureFields, &
          analysisThermoFields, &
          modelThermoFields, &
          sacaStateFields, sacaObsFields

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
     integer, allocatable, public :: nvert(:)                             ! number of vertical levels of each field

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
     procedure :: populate     => populate_subfields
     procedure :: delete       => delete_fields
     procedure :: read_file    => read_fields
     procedure :: write_file   => write_fields
     procedure :: serial_size  => serial_size
     procedure :: serialize    => serialize_fields
     procedure :: deserialize  => deserialize_fields
     procedure :: to_fieldset
     procedure :: from_fieldset
     procedure :: to_fieldset_ad
     !has
     generic, public :: has => has_field, has_fields
     procedure :: has_field
     procedure :: has_fields

     !get
     generic, public :: get => &
        get_data, &
        get_field_i1, get_field_i2, &
        get_field_r1, get_field_r2, &
        get_array_i1, get_array_i2, &
        get_array_r1, get_array_r2
     procedure :: &
        get_data, &
        get_field_i1, get_field_i2, &
        get_field_r1, get_field_r2, &
        get_array_i1, get_array_i2, &
        get_array_r1, get_array_r2

     !copy_to
     generic, public :: copy_to => &
        copy_to_other_fields_field, &
        copy_to_other_fields, &
        copy_to_other_pool_field, &
        copy_to_other_pool
     procedure :: &
        copy_to_other_fields_field, &
        copy_to_other_fields, &
        copy_to_other_pool_field, &
        copy_to_other_pool

     !copy_to_ad
     generic, public :: copy_to_ad => &
        copy_to_other_fields_field_ad, &
        copy_to_other_fields_ad, &
        copy_to_other_pool_field_ad, &
        copy_to_other_pool_ad
     procedure :: &
        copy_to_other_fields_field_ad, &
        copy_to_other_fields_ad, &
        copy_to_other_pool_field_ad, &
        copy_to_other_pool_ad

     !copy_from
     generic, public :: copy_from => &
        copy_from_other_fields_field, &
        copy_from_other_fields, &
        copy_from_other_pool_field, &
        copy_from_other_pool
     procedure :: &
        copy_from_other_fields_field, &
        copy_from_other_fields, &
        copy_from_other_pool_field, &
        copy_from_other_pool

     !push_back
     generic, public :: push_back => &
        push_back_other_fields_field, &
        push_back_other_fields, &
        push_back_other_pool_field, &
        push_back_other_pool
     procedure :: &
        push_back_other_fields_field, &
        push_back_other_fields, &
        push_back_other_pool_field, &
        push_back_other_pool


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
      "qc", "qi", "qr", "qs", "qg", "qh" ]
   character(len=MAXVARLEN) :: mpas_re_fields(3) = &
      [ character(len=MAXVARLEN) :: &
      "re_cloud", "re_ice  ", "re_snow " ]
   character(len=MAXVARLEN), parameter :: cellCenteredWindFields(2) = &
      [character(len=MAXVARLEN) :: &
       'uReconstructZonal', 'uReconstructMeridional']
   character(len=MAXVARLEN), parameter :: moistureFields(2) = &
      [character(len=MAXVARLEN) :: &
       'qv', 'spechum']
   character(len=MAXVARLEN), parameter :: analysisThermoFields(2) = &
      [character(len=MAXVARLEN) :: &
       'surface_pressure', 'temperature']
   character(len=MAXVARLEN), parameter :: modelThermoFields(4) = &
      [character(len=MAXVARLEN) :: &
       'qv', 'pressure', 'rho', 'theta']
   character(len=MAXVARLEN), parameter :: sacaStateFields(9) = &
      [character(len=MAXVARLEN) :: &
       'qv', 'qc', 'qi', 'qs', 'cldfrac', 'rho', 'temperature', 'pressure', 'xland']
   character(len=MAXVARLEN), parameter :: sacaObsFields(2) = &
      [character(len=MAXVARLEN) :: &
       'cldmask', 'brtemp']


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

    call self%populate()

    ! pre-determine number of vertical levels for each variables
    allocate(self % nvert(self % nf))
    do ivar = 1, self % nf
       self % nvert(ivar) = getVertLevels(self % subFields, self % fldnames(ivar))
    end do

    return

end subroutine create_fields

! ------------------------------------------------------------------------------

subroutine populate_subFields(self)

    implicit none
    class(mpas_fields), intent(inout) :: self

    call da_template_pool(self % geom, self % subFields, self % nf, self % fldnames)

end subroutine populate_subFields

! ------------------------------------------------------------------------------

subroutine delete_fields(self)

   implicit none
   class(mpas_fields), intent(inout) :: self
   integer :: ierr = 0

   if (allocated(self % fldnames)) deallocate(self % fldnames)
   if (allocated(self % fldnames_ci)) deallocate(self % fldnames_ci)
   if (allocated(self % nvert)) deallocate(self % nvert)

   call fckit_log%debug('--> delete_fields: deallocate subFields Pool')
   call delete_pool(self % subFields)

   call mpas_destroy_clock(self % clock, ierr)
   if ( ierr .ne. 0  ) then
      call fckit_log%info ('--> delete_fields deallocate clock failed')
   end if
   call fckit_log%debug('--> delete_fields done')

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

   call fckit_log%debug('--> copy_fields: copy subFields Pool')

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

   call fckit_log%debug('--> copy_fields done')

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
   logical :: Model2AnalysisVariableChange
   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer    :: pressure, pressure_base, pressure_p

   call fckit_log%debug('--> read_fields')
   call f_conf%get_or_die("date",str)
   sdate = str
   call datetime_set(sdate, vdate)

   call f_conf%get_or_die("filename",str)
   call swap_name_member(f_conf, str)
   temp_filename = str
   write(message,*) '--> read_fields: Reading ',trim(temp_filename)
   call fckit_log%debug(message)

   ! streamID (default: background)
   ! Name of the stream in streams.atmosphere or 'streams_file' associated with self%geom
   ! associated with this state.  Can be any string as long as it is included within the
   ! applicable streams.atmosphere file. Examples of stream names in the MPAS-JEDI distribution
   ! are 'background', 'analysis', 'ensemble', 'control', 'da_state'. Each of those streams has
   ! unique properties, including the MPAS fields that are read/written.
   streamID = 'background'
   if (f_conf%get("stream name", str)) then
     streamID = str
   end if

   ! temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   ierr = 0
   self % manager => self % geom % domain % streamManager
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(sdate,dateTimeString,1)
   write(message,*) '--> read_fields: dateTimeString: ',trim(dateTimeString)
   call fckit_log%debug(message)
   call mpas_set_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
   call mpas_set_clock_time(self % clock, local_time, MPAS_NOW)
   call mpas_set_clock_time(self % geom % domain % clock, local_time, MPAS_START_TIME)
   call mpas_expand_string(dateTimeString, -1, temp_filename, filename)
   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(message,*) '--> read_fields: Reading ',trim(filename)
   call fckit_log%debug(message)
   call MPAS_stream_mgr_read(self % manager, streamID=streamID, &
                           & when=dateTimeString, rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      write(message,*) '--> read_fields: MPAS_stream_mgr_read failed ierr=',ierr
      call abor1_ftn(message)
   end if

   ! Model2AnalysisVariableChange (default: true):
   ! indicates whether to transform from model fields (pressure_p, pressure_base, theta, qv)
   ! to analysis fields (temperature, specific_humidity).
   ! When streamID=='control' or 'saca_obs', the default value is changed to false.
   ! For example, the transform is not carried out when reading analysis fields directly (e.g.,
   ! background error standard deviation is read/written using streamID=='control') or
   ! when reading the obs-related fields for Non-Variational SAtellite-based Cloud Analysis (SACA).
   Model2AnalysisVariableChange = .True.
   if(streamID == 'control'.or. streamID == 'saca_obs') Model2AnalysisVariableChange = .False.
   if (f_conf%has("transform model to analysis")) then
      call f_conf%get_or_die("transform model to analysis", Model2AnalysisVariableChange)
   end if

   if(Model2AnalysisVariableChange) then
      !(1) diagnose pressure
      call mpas_pool_get_subpool(self % geom % domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_field(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_field(diag, 'pressure_base', pressure_base)
      call mpas_pool_get_field(diag, 'pressure', pressure)
      ngrid = self % geom % nCellsSolve
      pressure%array(:,1:ngrid) = pressure_base%array(:,1:ngrid) + pressure_p%array(:,1:ngrid)

      !(2) copy all to subFields & diagnose temperature
      call update_diagnostic_fields(self % geom % domain, self % subFields, self % geom % nCellsSolve)
   else
      call da_copy_all2sub_fields(self % geom % domain, self % subFields)
   endif

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
   write(message,*) '--> write_fields: ',trim(validitydate)
   call fckit_log%debug(message)
   call f_conf%get_or_die("filename",str)
   call swap_name_member(f_conf, str)
   temp_filename = str
   write(message,*) '--> write_fields: ',trim(temp_filename)
   call fckit_log%debug(message)
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
   write(message,*) 'check time --> write_fields: write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
   call fckit_log%debug(message)
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)

   self % manager => self % geom % domain % streamManager

   ! streamID (default: da_state)
   ! Name of the stream in streams.atmosphere or 'streams_file' associated with self%geom
   ! associated with this state.  Can be any string as long as it is included within the
   ! applicable streams.atmosphere file. Examples of stream names in the MPAS-JEDI distribution
   ! are 'background', 'analysis', 'ensemble', 'control', 'da_state'. Each of those streams has
   ! unique properties, including the MPAS fields that are read/written.
   streamID = 'da_state'
   if (f_conf%get("stream name", str)) then
     streamID = str
   end if

   call MPAS_stream_mgr_set_property(self % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)

   write(message,*) '--> write_fields: writing ',trim(filename)
   call fckit_log%debug(message)
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
   real(kind=RKIND), intent(out) :: pstat(3, nf)

   call da_gpnorm(self % subFields, self % geom % domain % dminfo, nf, pstat, fld_select = self % fldnames_ci(1:nf))

end subroutine gpnorm_

! ------------------------------------------------------------------------------

subroutine rms_(self, prms)

   implicit none
   class(mpas_fields),   intent(in)  :: self
   real(kind=RKIND), intent(out) :: prms

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
   real(kind=RKIND), intent(in)    :: zz

   call da_self_mult(self % subFields, zz)

end subroutine self_mult_

! ------------------------------------------------------------------------------

subroutine axpy_(self,zz,rhs)

   implicit none
   class(mpas_fields),   intent(inout) :: self
   real(kind=RKIND), intent(in)    :: zz
   class(mpas_fields),   intent(in)    :: rhs

   call da_axpy(self % subFields, rhs % subFields, zz, fld_select = self % fldnames_ci)

end subroutine axpy_

! ------------------------------------------------------------------------------

subroutine dot_prod_(self,fld,zprod)

   implicit none
   class(mpas_fields),    intent(in)    :: self, fld
   real(kind=RKIND),  intent(inout) :: zprod

   call da_dot_product(self % subFields, fld % subFields, self % geom % domain % dminfo, zprod)

end subroutine dot_prod_

!------------------------------------------------------------------------------

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
  real (kind=RKIND), dimension(:), pointer :: r1d_ptr
  real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr
  integer :: rhs_nCells, self_nCells, maxlevels, nlevels, jlev
  character(len=StrKIND) :: interp_type
  integer, allocatable :: rhsDims(:)

  interp_type = trim(rhs%geom%fields_to_fields_interp_type)

  if (interp_type == 'bump') then
    call initialize_bumpinterp(self%geom, rhs%geom, bumpinterp)
  else
    call initialize_uns_interp(self%geom, rhs%geom, unsinterp)
  endif

  ! Interpolate field from rhs mesh to self mesh using pre-calculated weights
  ! ----------------------------------------------------------------------------------
  maxlevels = rhs%geom%nVertLevelsP1
  rhs_nCells = rhs%geom%nCellsSolve
  self_nCells = self%geom%nCellsSolve

  allocate(interp_in(rhs_nCells, maxlevels))
  allocate(interp_out(self_nCells, maxlevels))

  call mpas_pool_begin_iteration(rhs%subFields)
  do while ( mpas_pool_get_next_member(rhs%subFields, poolItr) )
  ! TODO: For now, we skip the interpolation for integer type.
  if (poolItr % memberType == MPAS_POOL_FIELD .and. poolItr % dataType == MPAS_POOL_REAL ) then
    write(message,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
    call fckit_log%debug(message)
    ! Get the rhs fields out of the rhs%subFields pool and into an array that
    ! can be passed to the interpolator apply functions.
    ! (Two cases to cover: 1D/2D real variables)
    if (poolItr % nDims == 1) then
      nlevels = 1
      call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), r1d_ptr)
      interp_in(:,1) = r1d_ptr(1:rhs_nCells)
    else if (poolItr % nDims == 2) then
      rhsDims = getSolveDimSizes(rhs%subFields, poolItr%memberName)
      nlevels = rhsDims(1)
      if (nlevels > maxlevels) then
        write(message,*) '--> interpolate_fields: nlevels > maxlevels, ', nlevels, maxlevels
        call abor1_ftn(message)
      endif
      call mpas_pool_get_array(rhs%subFields, trim(poolItr % memberName), r2d_ptr)
      interp_in(1:rhs_nCells,1:nlevels) = transpose(r2d_ptr(1:nlevels,1:rhs_nCells))
    else
      write(message,*) '--> interpolate_fields: poolItr % nDims == ',poolItr % nDims,' not handled'
      call abor1_ftn(message)
    endif

    if (interp_type == 'bump') then
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
      call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r1d_ptr)
      r1d_ptr(1:self_nCells) = interp_out(:,1)
    else if (poolItr % nDims == 2) then
      call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), r2d_ptr)
      r2d_ptr(1:nlevels,1:self_nCells) = transpose (interp_out(1:self_nCells,1:nlevels))
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

   call bumpinterp%init(geom_from%f_comm,geom_from%afunctionspace,geom_to%afunctionspace,geom_from%nVertLevels)

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
   integer(c_size_t),intent(out) :: vsize !< Size

   ! Local variables
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: dimSizes(:)

   ! Initialize
   vsize = 0

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         dimSizes = getSolveDimSizes(self%subFields, poolItr%memberName)
         vsize = vsize + product(dimSizes)
         deallocate(dimSizes)
      endif
   enddo

end subroutine serial_size

! ------------------------------------------------------------------------------

subroutine serialize_fields(self, vsize, vect_inc)

   implicit none

   ! Passed variables
   class(mpas_fields),intent(in) :: self          !< Increment
   integer(c_size_t),intent(in) :: vsize          !< Size
   real(c_real_type),intent(out) :: vect_inc(vsize) !< Vector

   ! Local variables
   integer :: index, nvert, nhoriz, vv, hh
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: dimSizes(:)

   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
   integer, dimension(:), pointer :: i1d_ptr_a
   integer, dimension(:,:), pointer :: i2d_ptr_a

   ! Initialize
   index = 0

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         dimSizes = getSolveDimSizes(self%subFields, poolItr%memberName)
         if (poolItr % nDims == 1) then
            nhoriz = dimSizes(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i1d_ptr_a)
               do hh = 1, nhoriz
                  vect_inc(index + 1) = real(i1d_ptr_a(hh), kind=c_real_type)
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
            nvert = dimSizes(1)
            nhoriz = dimSizes(2)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     vect_inc(index + 1) = real(i2d_ptr_a(vv, hh), kind=c_real_type)
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
         deallocate(dimSizes)
      endif
   enddo

end subroutine serialize_fields

! --------------------------------------------------------------------------------------------------

subroutine deserialize_fields(self, vsize, vect_inc, index)

   implicit none

   ! Passed variables
   class(mpas_fields),intent(inout) :: self      !< Increment
   integer(c_size_t),intent(in) :: vsize         !< Size
   real(c_real_type),intent(in) :: vect_inc(vsize) !< Vector
   integer(c_size_t),intent(inout) :: index      !< Index

   ! Local variables
   integer :: nvert, nhoriz, vv, hh
   type (mpas_pool_iterator_type) :: poolItr
   integer, allocatable :: dimSizes(:)

   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
   integer, dimension(:), pointer :: i1d_ptr_a
   integer, dimension(:,:), pointer :: i2d_ptr_a

   call mpas_pool_begin_iteration(self%subFields)
   do while ( mpas_pool_get_next_member(self%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         dimSizes = getSolveDimSizes(self%subFields, poolItr%memberName)
         if (poolItr % nDims == 1) then
            nhoriz = dimSizes(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i1d_ptr_a)
               do hh = 1, nhoriz
                  i1d_ptr_a(hh) = nint ( vect_inc(index + 1) )
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
            nvert = dimSizes(1)
            nhoriz = dimSizes(2)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(self%subFields, trim(poolItr % memberName), i2d_ptr_a)
               do vv = 1, nvert
                  do hh = 1, nhoriz
                     i2d_ptr_a(vv, hh) = nint ( vect_inc(index + 1) )
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
         deallocate(dimSizes)
      endif
   enddo

end subroutine deserialize_fields

! has
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
   logical, allocatable :: has(:)
   allocate(has(size(fieldnames)))
   do i = 1, size(fieldnames)
      has(i) = self%has(fieldnames(i))
   end do
end function has_fields

! get
subroutine get_data(self, key, data)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   type(mpas_pool_data_type), pointer, intent(out) :: data
   if (self%has(key)) then
     data => pool_get_member(self % subFields, key, MPAS_POOL_FIELD)
   else
     write(message,*) 'self%get_data: field not present, ', key
     call abor1_ftn(message)
   end if
end subroutine get_data

subroutine get_field_i1(self, key, i1)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   type(field1DInteger), pointer, intent(out) :: i1
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   i1 => data%i1
end subroutine get_field_i1

subroutine get_field_i2(self, key, i2)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   type(field2DInteger), pointer, intent(out) :: i2
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   i2 => data%i2
end subroutine get_field_i2

subroutine get_field_r1(self, key, r1)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   type(field1DReal), pointer, intent(out) :: r1
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   r1 => data%r1
end subroutine get_field_r1

subroutine get_field_r2(self, key, r2)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   type(field2DReal), pointer, intent(out) :: r2
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   r2 => data%r2
end subroutine get_field_r2

subroutine get_array_i1(self, key, i1)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   integer, pointer, intent(out) :: i1(:)
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   i1 => data%i1%array
end subroutine get_array_i1

subroutine get_array_i2(self, key, i2)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   integer, pointer, intent(out) :: i2(:,:)
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   i2 => data%i2%array
end subroutine get_array_i2

subroutine get_array_r1(self, key, r1)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   real(kind=RKIND), pointer, intent(out) :: r1(:)
   type(mpas_pool_data_type), pointer :: data

   call self%get(key, data)
   r1 => data%r1%array
end subroutine get_array_r1

subroutine get_array_r2(self, key, r2)
   class(mpas_fields), intent(in) :: self
   character (len=*), intent(in) :: key
   real(kind=RKIND), pointer, intent(out) :: r2(:,:)
   type(mpas_pool_data_type), pointer :: data
   call self%get(key, data)
   r2 => data%r2%array
end subroutine get_array_r2


! all copy_to and copy_from methods eventually call
! copy_field_between_pools
subroutine copy_field_between_pools(from, fromKey, to, toKey)
  type(mpas_pool_type), pointer, intent(in) :: from
  type(mpas_pool_type), pointer, intent(inout) :: to
  character (len=*), intent(in) :: fromKey, toKey
  type(mpas_pool_data_type), pointer :: toData, fromData
  toData => pool_get_member(to, toKey, MPAS_POOL_FIELD)
  if (associated(toData)) then
    fromData => pool_get_member(from, fromKey, MPAS_POOL_FIELD)
    if (associated(fromData)) then
      if (associated(fromData%r1) .and. associated(toData%r1)) then
        toData%r1%array = fromData%r1%array
      else if (associated(fromData%r2) .and. associated(toData%r2)) then
        toData%r2%array = fromData%r2%array
      else if (associated(fromData%r3) .and. associated(toData%r3)) then
        toData%r3%array = fromData%r3%array
      else if (associated(fromData%i1) .and. associated(toData%i1)) then
        toData%i1%array = fromData%i1%array
      else if (associated(fromData%i2) .and. associated(toData%i2)) then
        toData%i2%array = fromData%i2%array
      else
        call abor1_ftn('copy_field_between_pools: data mismatch between to/from pools')
      end if
    else
      write(message,*) 'copy_field_between_pools: field not present in "from" pool, ', fromKey
      call abor1_ftn(message)
    end if
  else
    write(message,*) 'copy_field_between_pools: field not present in "to" pool, ', toKey
    call abor1_ftn(message)
  end if
end subroutine copy_field_between_pools

! copy_from
subroutine copy_from_other_pool_field(self, selfKey, otherPool, otherKey)
  class(mpas_fields), intent(inout) :: self
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  character (len=*), intent(in) :: selfKey, otherKey
  type(mpas_pool_data_type), pointer :: selfData, otherData
  call copy_field_between_pools(otherPool, otherKey, self%subFields, selfKey)
end subroutine copy_from_other_pool_field

subroutine copy_from_other_pool(self, key, otherPool)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  call self%copy_from(key, otherPool, key)
end subroutine copy_from_other_pool

subroutine copy_from_other_fields_field(self, selfKey, other, otherKey)
  class(mpas_fields), intent(inout) :: self
  class(mpas_fields), intent(in) :: other
  character (len=*), intent(in) :: selfKey, otherKey
  call self%copy_from(selfKey, other%subFields, otherKey)
end subroutine copy_from_other_fields_field

subroutine copy_from_other_fields(self, key, other)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  class(mpas_fields), intent(in) :: other
  call self%copy_from(key, other%subFields, key)
end subroutine copy_from_other_fields

! copy_to
subroutine copy_to_other_pool_field(self, selfKey, otherPool, otherKey)
  class(mpas_fields), intent(in) :: self
  type(mpas_pool_type), pointer, intent(inout) :: otherPool
  character (len=*), intent(in) :: selfKey, otherKey
  type(mpas_pool_data_type), pointer :: selfData, otherData
  call copy_field_between_pools(self%subFields, selfKey, otherPool, otherKey)
end subroutine copy_to_other_pool_field

subroutine copy_to_other_pool(self, key, otherPool)
  class(mpas_fields), intent(in) :: self
  character (len=*), intent(in) :: key
  type(mpas_pool_type), pointer, intent(inout) :: otherPool
  call self%copy_to(key, otherPool, key)
end subroutine copy_to_other_pool

subroutine copy_to_other_fields_field(self, selfKey, other, otherKey)
  class(mpas_fields), intent(in) :: self
  class(mpas_fields), intent(inout) :: other
  character (len=*), intent(in) :: selfKey, otherKey
  call self%copy_to(selfKey, other%subFields, otherKey)
end subroutine copy_to_other_fields_field

subroutine copy_to_other_fields(self, key, other)
  class(mpas_fields), intent(in) :: self
  character (len=*), intent(in) :: key
  class(mpas_fields), intent(inout) :: other
  call self%copy_to(key, other%subFields, key)
end subroutine copy_to_other_fields

! all copy_to_ad methods eventually call
! copy_field_between_pools_ad
subroutine copy_field_between_pools_ad(to, toKey, from, fromKey)
  type(mpas_pool_type), pointer, intent(inout) :: to
  type(mpas_pool_type), pointer, intent(in) :: from
  character (len=*), intent(in) :: fromKey, toKey
  type(mpas_pool_data_type), pointer :: toData, fromData
  toData => pool_get_member(to, toKey, MPAS_POOL_FIELD)
  if (associated(toData)) then
    fromData => pool_get_member(from, fromKey, MPAS_POOL_FIELD)
    if (associated(fromData)) then
      if (associated(fromData%r1) .and. associated(toData%r1)) then
        toData%r1%array = toData%r1%array + fromData%r1%array
      else if (associated(fromData%r2) .and. associated(toData%r2)) then
        toData%r2%array = toData%r2%array + fromData%r2%array
      else if (associated(fromData%r3) .and. associated(toData%r3)) then
        toData%r3%array = toData%r3%array + fromData%r3%array
      else
        call abor1_ftn('copy_field_between_pools_ad: data mismatch between to/from pools')
      end if
    else
      write(message,*) 'copy_field_between_pools_ad: field not present in "from" pool, ', fromKey
      call abor1_ftn(message)
    end if
  else
    write(message,*) 'copy_field_between_pools_ad: field not present in "to" pool, ', toKey
    call abor1_ftn(message)
  end if
end subroutine copy_field_between_pools_ad

! copy_to_ad
subroutine copy_to_other_pool_field_ad(self, selfKey, otherPool, otherKey)
  class(mpas_fields), intent(inout) :: self
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  character (len=*), intent(in) :: selfKey, otherKey
  type(mpas_pool_data_type), pointer :: selfData, otherData
  call copy_field_between_pools_ad(self%subFields, selfKey, otherPool, otherKey)
end subroutine copy_to_other_pool_field_ad

subroutine copy_to_other_pool_ad(self, key, otherPool)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  call self%copy_to_ad(key, otherPool, key)
end subroutine copy_to_other_pool_ad

subroutine copy_to_other_fields_field_ad(self, selfKey, other, otherKey)
  class(mpas_fields), intent(inout) :: self
  class(mpas_fields), intent(in) :: other
  character (len=*), intent(in) :: selfKey, otherKey
  call self%copy_to_ad(selfKey, other%subFields, otherKey)
end subroutine copy_to_other_fields_field_ad

subroutine copy_to_other_fields_ad(self, key, other)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  class(mpas_fields), intent(in) :: other
  call self%copy_to_ad(key, other%subFields, key)
end subroutine copy_to_other_fields_ad

! push_back
! all push_back methods eventually call
! pool_push_back_field_from_pool
subroutine pool_push_back_field_from_pool(to, toKey, from, fromKey)
  type(mpas_pool_type), pointer, intent(inout) :: to
  type(mpas_pool_type), pointer, intent(in) :: from
  character (len=*), intent(in) :: fromKey, toKey
  type(mpas_pool_data_type), pointer :: fromData
  type(field1DReal), pointer :: fieldr1
  type(field2DReal), pointer :: fieldr2
  type(field3DReal), pointer :: fieldr3
  type(field1DInteger), pointer :: fieldi1
  type(field2DInteger), pointer :: fieldi2
  fromData => pool_get_member(from, fromKey, MPAS_POOL_FIELD)
  if (associated(fromData)) then
    if (associated(fromData%r1)) then
      call mpas_duplicate_field(fromData%r1, fieldr1)
      fieldr1 % fieldName = toKey
      call mpas_pool_add_field(to, toKey, fieldr1)
    else if (associated(fromData%r2)) then
      call mpas_duplicate_field(fromData%r2, fieldr2)
      fieldr2 % fieldName = toKey
      call mpas_pool_add_field(to, toKey, fieldr2)
    else if (associated(fromData%r3)) then
      call mpas_duplicate_field(fromData%r3, fieldr3)
      fieldr3 % fieldName = toKey
      call mpas_pool_add_field(to, toKey, fieldr3)
    else if (associated(fromData%i1)) then
      call mpas_duplicate_field(fromData%i1, fieldi1)
      fieldi1 % fieldName = toKey
      call mpas_pool_add_field(to, toKey, fieldi1)
    else if (associated(fromData%i2)) then
      call mpas_duplicate_field(fromData%i2, fieldi2)
      fieldi2 % fieldName = toKey
      call mpas_pool_add_field(to, toKey, fieldi2)
    else
      call abor1_ftn('pool_push_back_field_from_pool: data type not supported')
    end if
  else
    write(message,*) 'pool_push_back_field_from_pool: field not present in "from" pool, ', fromKey
    call abor1_ftn(message)
  end if
end subroutine pool_push_back_field_from_pool

subroutine push_back_other_pool_field(self, selfKey, otherPool, otherKey)
  class(mpas_fields), intent(inout) :: self
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  character (len=*), intent(in) :: selfKey, otherKey
  type(mpas_pool_data_type), pointer :: selfData, otherData
  character(len=MAXVARLEN), allocatable :: fldnames(:)
  if (self%has(selfKey)) then
    write(message,*) 'push_back_other_pool_field: field already present in self, cannot push_back, ', selfKey
    call abor1_ftn(message)
  end if

  ! Add field to self%subFields pool
  call pool_push_back_field_from_pool(self%subFields, selfKey, otherPool, otherKey)

  ! Extend self%fldnames
  allocate(fldnames(self%nf+1))
  fldnames(1:self%nf) = self%fldnames(:)
  fldnames(self%nf+1) = trim(selfKey)
  self%nf = self%nf+1
  deallocate(self%fldnames)
  allocate(self%fldnames(self%nf))
  self%fldnames = fldnames
  deallocate(fldnames)
end subroutine push_back_other_pool_field

subroutine push_back_other_pool(self, key, otherPool)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  type(mpas_pool_type), pointer, intent(in) :: otherPool
  call self%push_back(key, otherPool, key)
end subroutine push_back_other_pool

subroutine push_back_other_fields_field(self, selfKey, other, otherKey)
  class(mpas_fields), intent(inout) :: self
  class(mpas_fields), intent(in) :: other
  character (len=*), intent(in) :: selfKey, otherKey
  call self%push_back(selfKey, other%subFields, otherKey)
end subroutine push_back_other_fields_field

subroutine push_back_other_fields(self, key, other)
  class(mpas_fields), intent(inout) :: self
  character (len=*), intent(in) :: key
  class(mpas_fields), intent(in) :: other
  call self%push_back(key, other%subFields, key)
end subroutine push_back_other_fields


subroutine to_fieldset(self, geom, vars, afieldset, include_halo, flip_vert_lev)

   implicit none

   class(mpas_fields),   intent(in)    :: self
   type(mpas_geom),      intent(in)    :: geom
   type(oops_variables), intent(in)    :: vars
   type(atlas_fieldset), intent(inout) :: afieldset
   logical,               intent(in)    :: include_halo
   logical,               intent(in)    :: flip_vert_lev

   integer :: jvar, nlevels
   real(kind=kind_real), pointer :: real_ptr(:,:)
   real(kind=RKIND), pointer :: r1d_ptr_a(:), r2d_ptr_a(:,:)
   integer, pointer :: i1d_ptr_a(:), i2d_ptr_a(:,:)
   logical :: var_found
   type(atlas_field) :: afield
   type(atlas_metadata) :: meta
   type(mpas_pool_iterator_type) :: poolItr

   type(mpas_pool_data_type), pointer :: data_aux
   type(field1DReal), pointer :: r1
   type(field2DReal), pointer :: r2
   type(field1DInteger), pointer :: i1
   type(field2DInteger), pointer :: i2
   integer :: nx, ilev, jlev
   integer :: j

   ! note:
   ! get-or-create atlas field, exch-halo in mpas-field, flip_vert_level,
   ! pass data/field to atlas, assign 'default, interp_type' in atlas's metadata

   if (include_halo) then
      nx=geom%nCells
   else
      nx=geom%nCellsSolve
   endif

   do jvar = 1,vars%nvars()
      var_found = .false.
      call mpas_pool_begin_iteration(self%subFields)
      do while (mpas_pool_get_next_member(self%subFields, poolItr))
         if (trim(vars%variable(jvar))==trim(poolItr%memberName)) then
            !
            ! Revealed a potetial bug in function getVertLevels(self%subFields, vars%variable(jvar))
            ! why getVertLevels .NE. 0 when nDims=1 ???
            ! Is there a purpose for this setup in getVertLevels ?
            !
            if (poolItr%nDims==1) then
               nlevels = 1
            else if (poolItr%nDims==2) then
               nlevels = getVertLevels(self%subFields, vars%variable(jvar))
            else if (poolItr%nDims==3) then
               call abor1_ftn('not implemented yet')
            end if

            if (afieldset%has_field(vars%variable(jvar))) then
               ! Get Atlas field
               afield = afieldset%field(vars%variable(jvar))
            else
               ! Create Atlas field
               if (include_halo) then
                  afield = geom%afunctionspace_incl_halo%create_field &
                       (name=vars%variable(jvar),kind=atlas_real(kind_real),levels=nlevels)
               else
                  afield = geom%afunctionspace%create_field &
                       (name=vars%variable(jvar),kind=atlas_real(kind_real),levels=nlevels)
               endif
               ! Add field
               call afieldset%add(afield)
            end if
            !
            call self%get_data(poolItr%memberName, data_aux)
            ! equiv. ref. data_aux = pool_get_member(self % subFields, poolItr%memberName, MPAS_POOL_FIELD)
            if (poolItr % dataType == MPAS_POOL_REAL) then
               if (poolItr%nDims==1) then
                  if (include_halo) then
                     r1 => data_aux%r1
                     call mpas_dmpar_exch_halo_field(r1)
                  endif
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r1d_ptr_a)
                  call afield%data(real_ptr)
                  real_ptr(1,1:nx) = real(r1d_ptr_a(1:nx), kind_real)
                  !
               else if (poolItr%nDims==2) then
                  if (include_halo) then
                     r2 => data_aux%r2
                     call mpas_dmpar_exch_halo_field(r2)
                  end if
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r2d_ptr_a)
                  call afield%data(real_ptr)
                  ! for CRTM: vertical level flip
                  do jlev = 1, nlevels
                     if (flip_vert_lev) then
                        ilev = nlevels - jlev + 1
                     else
                        ilev = jlev
                     endif
                     real_ptr(ilev,1:nx) = real(r2d_ptr_a(jlev,1:nx), kind_real)
                  enddo
               end if
            elseif (poolItr % dataType == MPAS_POOL_INTEGER) then
               if (poolItr%nDims==1) then
                  if (include_halo) then
                     i1 => data_aux%i1
                     call mpas_dmpar_exch_halo_field(i1)
                  end if
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), i1d_ptr_a)
                  call afield%data(real_ptr)
                  real_ptr(1,1:nx) = real(i1d_ptr_a(1:nx), kind_real)
               else if (poolItr%nDims==2) then
                  write(message,*) '--> fill_geovals: nDims == 2:  not handled for integers'
                  call abor1_ftn(message)
               end if
            else
               STOP 'poolItr % dataType neither real nor integer'
            endif


            meta = afield%metadata()
            if (poolItr % dataType == MPAS_POOL_REAL) then
               call meta%set('interp_type', 'default')
            elseif (poolItr % dataType == MPAS_POOL_INTEGER) then
               call meta%set('interp_type', 'integer')
            else
               call abor1_ftn('poolItr % dataType .NE. real OR integer, unexpected')
            endif

            ! Release pointer
            call afield%final()

            ! Set flag
            var_found = .true.
            exit
         end if   ! if (trim(vars%variable(jvar))==trim(poolItr%memberName))
      end do      ! mpas_pool_get_next_member
      if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
   end do         ! do jvar=1, vars%nvars()
end subroutine to_fieldset


subroutine from_fieldset(self, geom, vars, afieldset, include_halo, flip_vert_lev)

   implicit none

   class(mpas_fields),   intent(inout) :: self
   type(mpas_geom),      intent(in)    :: geom
   type(oops_variables), intent(in)    :: vars
   type(atlas_fieldset), intent(in)    :: afieldset
   logical,              intent(in)    :: include_halo
   logical,              intent(in)    :: flip_vert_lev

   integer :: jvar
   real(kind=kind_real), pointer :: real_ptr(:,:)
   real(kind=RKIND), pointer :: r1d_ptr_a(:), r2d_ptr_a(:,:)
   logical :: var_found
   type(atlas_field) :: afield
   type(mpas_pool_iterator_type) :: poolItr

   type(mpas_pool_data_type), pointer :: data_aux
   type(field1DReal), pointer :: r1
   type(field2DReal), pointer :: r2


   integer :: nlevels, nx, ilev, jlev
   ! Note:
   !   1. pass data from atlas_field to MPAS field
   !      a. real-type(afield) --> real/integer (MPAS type), bc atlas has only real
   !      b. flip back vertical levels
   !   2. no-halo exchange in this case for saber useage

   if (include_halo) then
      nx=geom%nCells
   else
      nx=geom%nCellsSolve
   endif

   do jvar = 1,vars%nvars()
      var_found = .false.
      call mpas_pool_begin_iteration(self%subFields)
      do while (mpas_pool_get_next_member(self%subFields, poolItr))
         if (trim(vars%variable(jvar))==trim(poolItr%memberName)) then
            ! Get atlas field
            afield = afieldset%field(vars%variable(jvar))
            ! Get MPAS data (pool_data_type)
            call self%get_data(poolItr%memberName, data_aux)
            nlevels = getVertLevels(self%subFields, vars%variable(jvar))
            if (poolItr % dataType == MPAS_POOL_REAL) then
               if (poolItr%nDims==1) then
                  call afield%data(real_ptr)
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r1d_ptr_a)
                  r1d_ptr_a(1:nx) = real(real_ptr(1,1:nx), RKIND)
               else if (poolItr%nDims==2) then
                  call afield%data(real_ptr)
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r2d_ptr_a)
                  ! for CRTM: vertical level flip
                  do jlev = 1, nlevels
                     if (flip_vert_lev) then
                        ilev = nlevels - jlev + 1
                     else
                        ilev = jlev
                     endif
                     r2d_ptr_a(jlev,1:nx) = real(real_ptr(ilev,1:nx), RKIND)
                  enddo
               end if
            elseif (poolItr % dataType == MPAS_POOL_INTEGER) then
                  write(message,*) 'from_fieldset error: integers are not handled here'
                  call abor1_ftn(message)
            else
               call abor1_ftn('poolItr % dataType neither real nor integer')
            endif

            ! Release pointer
            call afield%final()

            ! Set flag
            var_found = .true.
            exit
         end if
      end do
      if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
   end do
end subroutine from_fieldset

! ------------------------------------------------------------------------------

subroutine to_fieldset_ad (self, geom, vars, afieldset, include_halo, flip_vert_lev)

   implicit none

   class(mpas_fields),   intent(inout) :: self
   type(mpas_geom),      intent(in)    :: geom
   type(oops_variables), intent(in)    :: vars
   type(atlas_fieldset), intent(in)    :: afieldset
   logical,              intent(in)    :: include_halo
   logical,              intent(in)    :: flip_vert_lev

   integer :: jvar
   real(kind=kind_real), pointer :: real_ptr(:,:)
   real(kind=RKIND), pointer :: r1d_ptr_a(:), r2d_ptr_a(:,:)
   logical :: var_found
   type(atlas_field) :: afield
   type(mpas_pool_iterator_type) :: poolItr

   type(mpas_pool_data_type), pointer :: data_aux
   type(field1DReal), pointer :: r1
   type(field2DReal), pointer :: r2

   integer :: nlevels, nx, ilev, jlev
   ! Note:
   !   this is the adjoint of to_fieldset, so reverse the code instruction order
   !   1. do assignment (Y_atlas=X_mpas --> X_mpas* = X* + Y* ) in vertical level flip
   !   2. exch_halo_ad :  AD code for MPI exch_halo

   if (include_halo) then
      nx=geom%nCells
   else
      nx=geom%nCellsSolve
      call abor1_ftn('Error in subroutine to_fieldset_ad, must have halo here')
   endif

   do jvar = 1,vars%nvars()
      var_found = .false.
      call mpas_pool_begin_iteration(self%subFields)
      do while (mpas_pool_get_next_member(self%subFields, poolItr))
         if (trim(vars%variable(jvar))==trim(poolItr%memberName)) then
            ! Get atlas field
            afield = afieldset%field(vars%variable(jvar))
            ! Get MPAS data (pool_data_type)
            call self%get_data(poolItr%memberName, data_aux)
            nlevels = getVertLevels(self%subFields, vars%variable(jvar))
            if (poolItr % dataType == MPAS_POOL_REAL) then
               if (poolItr%nDims==1) then
                  call afield%data(real_ptr)
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r1d_ptr_a)
                  r1d_ptr_a(1:nx) = r1d_ptr_a(1:nx) + real(real_ptr(1,1:nx), RKIND)
                  r1 => data_aux%r1
                  call mpas_dmpar_exch_halo_adj_field(r1)
               else if (poolItr%nDims==2) then
                  call afield%data(real_ptr)
                  call mpas_pool_get_array(self%subFields, trim(poolItr%memberName), r2d_ptr_a)
                  ! for CRTM: vertical level flip
                  do jlev = nlevels, 1, -1
                     if (flip_vert_lev) then
                        ilev = nlevels - jlev + 1
                     else
                        ilev = jlev
                     endif
                     r2d_ptr_a(jlev,1:nx) = r2d_ptr_a(jlev,1:nx) + real(real_ptr(ilev,1:nx), RKIND)
                  enddo
                  r2 => data_aux%r2
                  call mpas_dmpar_exch_halo_adj_field(r2)
               end if
            elseif (poolItr % dataType == MPAS_POOL_INTEGER) then
               call abor1_ftn('to_fieldset_ad integer should not happen, stop')
            endif
            ! Release pointer
            call afield%final()

            ! Set flag
            var_found = .true.
            exit
         end if
      end do
      if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
   end do
end subroutine to_fieldset_ad


end module mpas_fields_mod
