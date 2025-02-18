! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_model_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use duration_mod
use mpas_geom_mod
use mpas_fields_mod
use mpas_trajectories

use mpas_derived_types
use mpas_framework
use mpas_kind_types
use mpas_subdriver
use atm_core
use mpas_stream_manager
use mpas4da_mod
use mpas_constants, only : rgas, cp

use kinds, only : kind_real

implicit none
private
public :: mpas_model, & 
        & model_setup, model_delete, &
        & model_prepare_integration, model_prepare_integration_tl, model_prepare_integration_ad, &
        & model_propagate, model_propagate_tl, model_propagate_ad, &
        & model_prop_traj, model_wipe_traj, &
        & mpas_model_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: mpas_model
   ! GD: for now, model is just using the full geom structure
   ! we update the fields for time integration using subfield
   ! from mpas_fields
   type (domain_type), pointer :: domain 
   type (core_type), pointer :: corelist
   real (kind=RKIND) :: dt
end type mpas_model

integer, parameter    :: max_string=8000
character(max_string) :: message

#define LISTED_TYPE mpas_model

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include <oops/util/linkedList_c.f>

! ------------------------------------------------------------------------------

subroutine model_setup(self, geom, f_conf)

!   use fckit_mpi_module, only: fckit_mpi_comm

   implicit none
   type(fckit_configuration), intent(in)    :: f_conf !< fckit config
   type(mpas_model),          intent(inout) :: self  ! should I put intent on these?
!   type(mpas_model), target                 :: model  ! should I put intent on these?
   type(mpas_geom)        :: geom

   character(len=:), allocatable :: str
   character(len=20) :: ststep
   type(duration) :: dtstep

   real (kind=RKIND), pointer :: config_dt
   character (len=StrKIND), pointer :: config_start_time
   character (len=StrKIND), pointer :: config_restart_timestamp_name
   character (len=StrKIND), pointer :: config_run_duration
   character (len=StrKIND), pointer :: config_stop_time
   character (len=StrKIND) :: startTimeStamp

!   type(fckit_mpi_comm) :: f_comm
   
   call fckit_log%info('===> model_setup')
#define ModelMPAS_setup
#ifdef ModelMPAS_setup
   self % corelist => geom % corelist
   self % domain => geom % domain
!   f_comm = fckit_mpi_comm()
!   !> MPAS subdriver
!   call mpas_init( self % corelist, self % domain, external_comm=f_comm%communicator() )
!
!   if (associated(self % domain)) then
!       call fckit_log%debug('inside model: model % domain associated')
!   end if
!   if (associated(self % corelist)) then
!       call fckit_log%debug('inside model: model % corelist associated')
!   else
!       call fckit_log%debug('inside model: model % corelist not associated')
!   end if

   ! GD:  we need to update some parameters here regarding the yaml namelist file of oops.
   ! Also, we can add new DA parameters in the MPAS configs file if needed.
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_dt', config_dt)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_start_time', config_start_time)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_run_duration', config_run_duration)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_stop_time', config_stop_time)
   write(message,*) 'config_dt: ',config_dt
   call fckit_log%debug(message)
   write(message,*) 'config_start_time: ',trim(config_start_time)
   call fckit_log%debug(message)
   write(message,*) 'config_restart_timestamp_name: ',trim(config_restart_timestamp_name)
   call fckit_log%debug(message)
   write(message,*) 'config_run_duration: ',trim(config_run_duration)
   call fckit_log%debug(message)
   write(message,*) 'config_stop_time: ',trim(config_stop_time)
   call fckit_log%debug(message)

   ! GD: needs a converter from oops yaml file format to mpas if the yaml file drives MPAS
   ! otherwise mpas namelist file can be used.
   call f_conf%get_or_die("tstep",str)
   ststep = str
   dtstep = trim(ststep)
   write(message,*) 'ststep: ', ststep
   call fckit_log%debug(message)
   self % dt = config_dt !real(duration_seconds(dtstep),kind_real)
   ! call f_conf%get_or_die("dstep",str)
   ! dstep = str
   ! set clock here MPAS_STARTING, MPAS_NOW, MPAS_STOP_TIME
   ! config_start_time = '2010-10-23_00:00:00'
   ! config_run_duration = '0_02:00:00'
#endif

end subroutine model_setup

! ------------------------------------------------------------------------------

subroutine model_delete(self)

   implicit none
   type(mpas_model) :: self

   ! For now, all the structure is hold by geom
   call fckit_log%info('===> model_delete')
!   if ((associated(self % corelist)).and.(associated(self % domain))) then
!      call fckit_log%debug('===> delete model corelist and domain')
!      call mpas_timer_set_context( self % domain )
!      call mpas_finalize(self % corelist, self % domain)
!   end if

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_prepare_integration(self, jedi_state)

   implicit none

   type(mpas_model) :: self
   type(mpas_fields) :: jedi_state
   logical, pointer :: config_do_restart, config_do_DAcycling, config_dt
   integer :: ierr = 0
   real (kind=RKIND), pointer :: dt
   type (block_type), pointer :: block

   character(len=StrKIND) :: startTimeStamp, stopTimeStamp, nowTimeStamp

   type (mpas_pool_type), pointer :: state
   type (mpas_pool_type), pointer :: mesh
   type (mpas_pool_type), pointer :: diag
   type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
   type (field2DReal), pointer :: uReconstructZonal, uReconstructMeridional
   character (len=StrKIND), pointer :: xtime
   character (len=StrKIND), pointer :: config_run_duration
   type (MPAS_Time_Type) :: startTime, stopTime, nowTime
   type (MPAS_Timeinterval_Type) ::  runDuration
   type (MPAS_Alarm_type), pointer :: alarmPtr

   call fckit_log%info('===> model_prepare_integration')

   !--------------------------
   ! GD: the present design relies on the hypothesis that we run
   ! the model member sequentially using the geom structure
   ! In the reverse case, we will need to create locally domain and 
   ! corelist by calling first the mpas_subdriver for example 
   ! like it is done in geom
   !----------------------------------------------------------------

   ! here or where increment is computed
   call da_copy_sub2all_fields(jedi_state % geom, jedi_state % subFields)

!#ifdef odelMPAS_prepare
   !-------------------------------------------------------------------
   ! WIND processing U and V resconstruct to the edges
   ! not parallel yet, routine initially from DART
   !-------------------------------------------------------------------

   ! here or where increment is computed
   !call mpas_pool_get_field(jedi_state % subfields, 'eastward_wind', uReconstructZonal)
   !call mpas_pool_get_field(jedi_state % subfields, 'northward_wind', uReconstructMeridional)
   !call mpas_pool_get_field(state, 'u', u_field, 1)
   !call uv_cell_to_edges(self%domain, uReconstructZonal, uReconstructMeridional, u_field, & 
   !                  & jedi_state%geom%lonCell, jedi_state%geom%latCell, &
   !                  & jedi_state%geom%nCellsGlobal, jedi_state%geom%edgeNormalVectors, &
   !                  & jedi_state%geom%nEdgesOnCell, jedi_state%geom%edgesOnCell, jedi_state%geom%nVertLevels)

   !-------------------------------------------------------------------
   ! update domain % clock using mpas_fields clock and config files
   !-------------------------------------------------------------------
   startTime = mpas_get_clock_time(jedi_state % clock, MPAS_START_TIME, ierr)
   call mpas_set_clock_time(self % domain % clock, startTime, MPAS_START_TIME)
   call mpas_get_time(startTime, dateTimeString=startTimeStamp) ! needed by xtime later
   call mpas_pool_get_config(self % domain % blocklist % configs,'config_run_duration',config_run_duration)
   call mpas_set_timeInterval(runDuration, timeString=config_run_duration,ierr=ierr)
   stopTime = startTime + runDuration
   call mpas_set_clock_time(self % domain % clock, stopTime, MPAS_STOP_TIME)
   call mpas_get_time(stopTime, dateTimeString=stopTimeStamp)
   write(message,*) 'MPAS_START_TIME, MPAS_STOP_TIME: ',trim(startTimeStamp),trim(stopTimeStamp)
   call fckit_log%debug(message)
 !--
   nowTime = mpas_get_clock_time(jedi_state % clock, MPAS_NOW, ierr)
   call mpas_get_time(nowTime, dateTimeString=nowTimeStamp)
   write(message,*) 'MPAS_NOW from jedi_state % clock: ',trim(nowTimeStamp)
   call fckit_log%debug(message)
   nowTime = mpas_get_clock_time(self % domain % clock, MPAS_NOW, ierr)
   call mpas_get_time(nowTime, dateTimeString=nowTimeStamp)
   write(message,*) 'MPAS_NOW from self % domain % clock: ',trim(nowTimeStamp)
   call fckit_log%debug(message)
!-- set xtime as startTimeStamp
   call mpas_pool_get_array(self % domain % blocklist % allFields, 'xtime', xtime, 1)
   write(message,*) 'xtime_old=',xtime
   call fckit_log%debug(message)
   xtime = startTimeStamp
   write(message,*) 'xtime_new=',xtime
   call fckit_log%debug(message)

   !--------------------------------------------------------------------
   ! Computation of theta_m from theta and qv
   ! Computation of rho_zz from rho / zz from updated rho
   ! Recoupling of all the variables
   !--------------------------------------------------------------------
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_do_restart', config_do_restart)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_do_DAcycling', config_do_DAcycling)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_dt', dt)
   config_do_restart = .True.
   config_do_DAcycling = .True.

!--- TODO: clean-up this subroutine !!
!   call mpas_pool_get_subpool(self % domain % blocklist % structs,'state',state)
!   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'diag', diag)
!   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'mesh', mesh)
!   call atm_compute_output_diagnostics(state, 1, diag, mesh)

!   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
!   call mpas_pool_get_field(state, 'u', u_field, 1)
!   call mpas_dmpar_exch_halo_field(u_field)
!   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'diag', diag)
!   call mpas_pool_get_field(diag, 'rho', u_field, 1)
!   call mpas_dmpar_exch_halo_field(u_field)
!   call mpas_pool_get_field(diag, 'theta', u_field, 1)
!   call mpas_dmpar_exch_halo_field(u_field)

   block => self % domain % blocklist
#define ModelMPAS_prepare
#ifdef ModelMPAS_prepare
   !! Below: "clock" pointer comes from atm_core
   !! "clock" needs to be associated with self % domain % clock in order to initialize alarms correctly
   !! in atm_mpas_init_block
   clock => self % domain % clock 

   !Destroy alarms that were previously initialized (will be re-initialized in atm_mpas_init_block)
   alarmPtr => clock % alarmListHead
   do while (associated(alarmPtr))
      clock % alarmListHead => alarmPtr % next
      deallocate(alarmPtr)
      alarmPtr => clock % alarmListHead
   end do

   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)

      ! GD: if we do cycling in atm_mpas_init_block we propably need to avoid recomputing wind to the 
      ! mass center. (avoiding doing twice ... adding a flag). OK for now, probably same problem with dart 
      call atm_mpas_init_block(self % domain % dminfo, self % domain % streamManager, block, mesh, self % dt)
      call mpas_pool_get_array(state, 'xtime', xtime, 1)
      xtime = startTimeStamp
      block => block % next
   end do

   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'diag', diag)
   call mpas_pool_get_field(diag, 'pv_edge', pv_edge_field)
   call mpas_dmpar_exch_halo_field(pv_edge_field)
   call mpas_pool_get_field(diag, 'ru', ru_field)
   call mpas_dmpar_exch_halo_field(ru_field)
   call mpas_pool_get_field(diag, 'rw', rw_field)
   call mpas_dmpar_exch_halo_field(rw_field)
#endif

end subroutine model_prepare_integration

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_ad(self, inc)

   implicit none
   type(mpas_model)  :: self
   type(mpas_fields) :: inc

   call fckit_log%info ('===> model_prepare_integration_ad')

end subroutine model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_tl(self, inc)

   implicit none
   type(mpas_model)  :: self
   type(mpas_fields) :: inc

   call fckit_log%info ('===> model_prepare_integration_tl')

end subroutine model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine model_propagate(self, jedi_state)

   implicit none
   type(mpas_model)  :: self
   type(mpas_fields) :: jedi_state

   type (mpas_pool_type), pointer :: state, diag, mesh
   real (kind=kind_real) :: dt
   integer :: itimestep

   call fckit_log%info ('===> model_propagate')

   ! Get state, diag, and mesh from domain
   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'diag', diag)
   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'mesh', mesh)

   itimestep = 1
   ! Perform single time step, including
   !  update of mpas wind from the edges to center
   call atm_do_timestep(self % domain, self % dt, itimestep)
   call mpas_pool_shift_time_levels(state)
   call mpas_advance_clock(self % domain % clock)
   call mpas_advance_clock(jedi_state % clock) !Advance jedi_state % clock, Too

   ! TODO: GD: can be here or needs to go probably somewhere else as a postprocessing of a forecast
   ! update theta et rho: postprocess is implemented at the C++ level now and needs to be
   ! interfaced with fortran.
   !if ( mpas_is_clock_time(self % domain % clock) ) then 
      !(1) diagnose theta, rho, pressure
      call atm_compute_output_diagnostics(state, 1, diag, mesh)

      !(2) copy all to subFields & diagnose temperature
      call update_diagnostic_fields(jedi_state % geom, jedi_state % subFields, jedi_state % geom % nCellsSolve)
   !end if

end subroutine model_propagate

! ------------------------------------------------------------------------------

subroutine model_propagate_ad(self, inc, traj)
   implicit none
   type(mpas_model)      :: self
   type(mpas_fields)     :: inc
   type(mpas_trajectory) :: traj
   call fckit_log%info ('===> model_propagate_ad')
end subroutine model_propagate_ad

! ------------------------------------------------------------------------------

subroutine model_propagate_tl(self, inc, traj)
   implicit none
   type(mpas_model)      :: self
   type(mpas_fields)     :: inc
   type(mpas_trajectory) :: traj
   call fckit_log%info ('===> model_propagate_tl')
end subroutine model_propagate_tl

! ------------------------------------------------------------------------------

subroutine model_prop_traj(self, jedi_state, traj)

   implicit none
   type(mpas_model)      :: self
   type(mpas_fields)     :: jedi_state
   type(mpas_trajectory) :: traj

   call fckit_log%info ('===> model_prop_traj in mpas_model_mod.F90')
   call set_traj(traj,jedi_state)

end subroutine model_prop_traj

! ------------------------------------------------------------------------------

subroutine model_wipe_traj(traj)

   implicit none
   type(mpas_trajectory) :: traj

   call fckit_log%info ('===> model_wipe_traj')

end subroutine model_wipe_traj

! ------------------------------------------------------------------------------

end module mpas_model_mod
