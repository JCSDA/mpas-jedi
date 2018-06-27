! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_model_mod

use iso_c_binding
use config_mod
use duration_mod
use mpas_geom_mod
use mpas_fields_mod
use mpas_trajectories

use mpas_derived_types
use mpas_framework
use mpas_kind_types
use mpas_subdriver
use atm_core
use mpi 
use mpas_stream_manager
use mpas4da_mod
use mpas_kinds, only : kind_real

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
   ! from mpas_field
   type (domain_type), pointer :: domain 
   type (core_type), pointer :: corelist
   real (kind=kind_real) :: dt
end type mpas_model

#define LISTED_TYPE mpas_model

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine model_setup(self, geom, c_conf)

   implicit none
   type(c_ptr), intent(in) :: c_conf !< pointer to object of class Config
   !type(mpas_model), target :: model  ! should I put intent on these?
   type(mpas_model), intent(inout) :: self  ! should I put intent on these?
   type(mpas_geom)        :: geom

   character(len=20) :: ststep
   type(duration) :: dtstep

   real (kind=kind_real), pointer :: config_dt
   character (len=StrKIND), pointer :: config_start_time
   character (len=StrKIND), pointer :: config_restart_timestamp_name
   character (len=StrKIND), pointer :: config_run_duration
   character (len=StrKIND), pointer :: config_stop_time
   character (len=StrKIND) :: startTimeStamp
   
   write(*,*) "---- Inside of Sub. model_setup ----"
#define ModelMPAS_setup
#ifdef ModelMPAS_setup
   self % corelist => geom % corelist
   self % domain => geom % domain

   ! GD:  we need to update some parameters here regarding the json namelist file of oops.
   ! Also, we can add new DA parameters in the MPAS configs file if needed.
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_dt', config_dt)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_start_time', config_start_time)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_run_duration', config_run_duration)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_stop_time', config_stop_time)
   write(*,*)'config_dt: ',config_dt
   write(*,*)'config_start_time: ',trim(config_start_time)
   write(*,*)'config_restart_timestamp_name: ',trim(config_restart_timestamp_name)
   write(*,*)'config_run_duration: ',trim(config_run_duration)
   write(*,*)'config_stop_time: ',trim(config_stop_time)
   write(0,*)'geom % nCellsGlobal: ',geom % nCellsGlobal
   write(0,*)'geom % nCells: ',geom % nCells
   write(0,*)'geom % nCellsSolve: ',geom % nCellsSolve



   ! GD: needs a converter from oops json file format to mpas if the json file drives MPAS
   ! otherwise mpas namelist file can be used.
   ststep = config_get_string(c_conf,len(ststep),"tstep")
   dtstep = trim(ststep)
   write(0,*)'ststep: ', ststep
   self % dt = config_dt !real(duration_seconds(dtstep),kind_real)
   ! dstep = config_get_string(c_conf,len(dtstep),"dstep")
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
   write(*,*)'===> model_delete'

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_prepare_integration(self, flds)

   implicit none

   type(mpas_model) :: self
   type(mpas_field) :: flds
   logical, pointer :: config_do_restart, config_do_DAcycling, config_dt
   integer :: ierr = 0
   real (kind=kind_real), pointer :: dt
   type (block_type), pointer :: block

   character(len=StrKIND) :: startTimeStamp, stopTimeStamp

   type (mpas_pool_type), pointer :: state
   type (mpas_pool_type), pointer :: mesh
   type (mpas_pool_type), pointer :: diag
   type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
   type (field2DReal), pointer :: uReconstructZonal, uReconstructMeridional
   character (len=StrKIND), pointer :: xtime
   character (len=StrKIND), pointer :: config_run_duration
   type (MPAS_Time_Type) :: startTime, stopTime
   type (MPAS_Timeinterval_Type) ::  runDuration

   write(*,*)'===> model_prepare_integration'

   !————————————————————————————————————————————---------------------------
   ! GD: the present design relies on the hypothesis that we run
   ! the model member sequentially using the geom structure
   ! In the reverse case, we will need to create locally domain and 
   ! corelist by calling first the mpas_subdriver for example 
   ! like it is done in geom
   !————————————————————————————————————————————---------------------------

   ! here or where increment is computed
    call da_copy_sub2all_fields(self % domain, flds % subFields)   
    call da_copy_sub2all_fields(self % domain, flds % auxFields)   

!#define ModelMPAS_prepare
#ifdef ModelMPAS_prepare
   !-------------------------------------------------------------------
   ! WIND processing U and V resconstruct to the edges
   ! not parallel yet, routine initially from DART
   !-------------------------------------------------------------------

   ! here or where increment is computed
   !call mpas_pool_get_field(flds % subfields, 'uReconstructZonal', uReconstructZonal)
   !call mpas_pool_get_field(flds % subfields, 'uReconstructMeridional', uReconstructMeridional)
   !call mpas_pool_get_field(state, 'u', u_field, 1)
   !call uv_cell_to_edges(self%domain, uReconstructZonal, uReconstructMeridional, u_field, & 
   !                  & flds%geom%lonCell, flds%geom%latCell, &
   !                  & flds%geom%nCellsGlobal, flds%geom%edgeNormalVectors, &
   !                  & flds%geom%nEdgesOnCell, flds%geom%edgesOnCell, flds%geom%nVertLevels)

   !——————————————————————————————————----------------------------------
   ! update domain % clock using mpas_field clock and config files
   !——————————————————————————————————---------------------------------
   startTime = mpas_get_clock_time(flds % clock, MPAS_START_TIME, ierr)
   call mpas_set_clock_time(self % domain % clock, startTime, MPAS_START_TIME)
   call mpas_get_time(startTime, dateTimeString=startTimeStamp) ! needed by xtime later
   call mpas_pool_get_config(self % domain % blocklist % configs,'config_run_duration',config_run_duration)
   call mpas_set_timeInterval(runDuration, timeString=config_run_duration,ierr=ierr)
   stopTime = startTime + runDuration
   call mpas_set_clock_time(self % domain % clock, stopTime, MPAS_STOP_TIME)
   call mpas_get_time(startTime, dateTimeString=stopTimeStamp)
   write(*,*)'MPAS_START_TIME, MPAS_STOP_TIME: ',trim(startTimeStamp),trim(stopTimeStamp)

   !——————————————————————————————————----------------------------------
   ! Computation of theta_m from theta and qv
   ! Computation of rho_zz from rho / zz from updated rho
   ! Recoupling of all the variables
   !——————————————————————————————————----------------------------------
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_do_restart', config_do_restart)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_do_DAcycling', config_do_DAcycling)
   call mpas_pool_get_config(self % domain % blocklist % configs, 'config_dt', dt)
   config_do_restart = .True.
   config_do_DAcycling = .True.

   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
   !call mpas_pool_get_field(state, 'u', u_field, 1)
   !call mpas_dmpar_exch_halo_field(u_field)

   block => self % domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
write(*,*) 'step 8'
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

subroutine model_prepare_integration_ad(self, flds)

   implicit none
   type(mpas_model) :: self
   type(mpas_field) :: flds

   write(*,*)'===> model_prepare_integration_ad'

end subroutine model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_tl(self, flds)

   implicit none
   type(mpas_model) :: self
   type(mpas_field) :: flds

   write(*,*)'===> model_prepare_integration_tl'

end subroutine model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine model_propagate(self, flds)

   implicit none
   type(mpas_model) :: self
   type(mpas_field) :: flds

   type (mpas_pool_type), pointer :: state, diag, mesh
   real (kind=kind_real) :: dt
   integer :: itimestep
!--bjj clock test
   character (len=StrKIND) :: tmpTimeStamp
   type (MPAS_Time_Type) :: tmpTime
   integer :: ierr = 0

   write(*,*)'===> model_propagate'

   itimestep = 1
   call atm_do_timestep(self % domain, self % dt, itimestep)
   call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
   call mpas_pool_shift_time_levels(state)
   call mpas_advance_clock(self % domain % clock)
   call mpas_advance_clock(flds % clock) !Advance flds % clock, Too
!--bjj clock test
   tmpTime = mpas_get_clock_time(self % domain % clock, MPAS_NOW, ierr)
   call mpas_get_time(tmpTime, dateTimeString=tmpTimeStamp) 
   write(*,*) '== BJJ in sub model_propagate time= ',trim(tmpTimeStamp)

   ! TODO: GD: can be here or needs to go probably somewhere else as a postprocessing of a forecast
   ! update theta et rho: postprocess is implemented at the C++ level now and needs to be
   ! interfaced with fortran.
   !if ( mpas_is_clock_time(self % domain % clock) ) then 
      !call mpas_pool_get_subpool(self % domain % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(self % domain % blocklist % structs, 'diag', diag)
      call mpas_pool_get_subpool(self % domain % blocklist % structs, 'mesh', mesh)
      call atm_compute_restart_diagnostics(state, 1, diag, mesh)
      call atm_compute_output_diagnostics(state, 1, diag, mesh) !--> add pressure
      ! update mpas wind from the edges to center
      ! copy from allfields to subfields
      call da_copy_all2sub_fields(self % domain, flds % subFields) 
      call da_copy_all2sub_fields(self % domain, flds % auxFields) 
   !end if

end subroutine model_propagate

! ------------------------------------------------------------------------------

subroutine model_propagate_ad(self, flds, traj)

   implicit none

   type(mpas_model)      :: self
   type(mpas_field)      :: flds
   type(mpas_trajectory) :: traj

   write(*,*)'===> model_prepare_integration_ad'

end subroutine model_propagate_ad

! ------------------------------------------------------------------------------

subroutine model_propagate_tl(self, flds, traj)

   implicit none
   type(mpas_model)      :: self
   type(mpas_field)      :: flds
   type(mpas_trajectory) :: traj

   write(*,*)'===> model_propagate_tl'

end subroutine model_propagate_tl

! ------------------------------------------------------------------------------

subroutine model_prop_traj(self, flds, traj)

   implicit none
   type(mpas_model)      :: self
   type(mpas_field)      :: flds
   type(mpas_trajectory) :: traj

   write(*,*)'===> model_prop_traj in mpas_model_mod.F90'
   call set_traj(traj,flds)

end subroutine model_prop_traj

! ------------------------------------------------------------------------------

subroutine model_wipe_traj(traj)

   implicit none
   type(mpas_trajectory) :: traj

   write(*,*)'===> model_wipe_traj'

end subroutine model_wipe_traj

! ------------------------------------------------------------------------------

end module mpas_model_mod
