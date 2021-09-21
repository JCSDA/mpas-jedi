! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_mpas_model_setup(c_conf, c_key_geom, c_key_self) bind (c,name='mpas_model_setup_f90')

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use duration_mod
use mpas_model_mod
use mpas_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr), value, intent(in) :: c_conf     !< pointer to object of class Config

type(mpas_model), pointer :: model
type(mpas_geom), pointer :: geom
type(fckit_configuration) :: f_conf

call mpas_geom_registry%get(c_key_geom, geom)
call mpas_model_registry%init()
call mpas_model_registry%add(c_key_self)
call mpas_model_registry%get(c_key_self, model)

f_conf = fckit_configuration(c_conf)
call model_setup(model, geom, f_conf)

end subroutine c_mpas_model_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_model_delete(c_key_self) bind (c,name='mpas_model_delete_f90')

use mpas_model_mod
use iso_c_binding

implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_model), pointer :: self

call mpas_model_registry%get(c_key_self, self)
call model_delete(self)
call mpas_model_registry%remove(c_key_self)

end subroutine c_mpas_model_delete

! ------------------------------------------------------------------------------

subroutine c_mpas_model_prepare_integration(c_key_self, c_key_state) &
         & bind(c,name='mpas_model_prepare_integration_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model fields

type(mpas_model), pointer :: self
type(mpas_fields), pointer :: state

call mpas_fields_registry%get(c_key_state,state)
call mpas_model_registry%get(c_key_self, self)

call model_prepare_integration(self, state)

end subroutine c_mpas_model_prepare_integration

! ------------------------------------------------------------------------------

subroutine c_mpas_model_prepare_integration_ad(c_key_self, c_key_incr) &
           bind(c,name='mpas_model_prepare_integration_ad_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model fields

type(mpas_model), pointer :: self
type(mpas_fields), pointer :: inc

call mpas_model_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_incr,inc)

call model_prepare_integration_ad(self, inc)

end subroutine c_mpas_model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine c_mpas_model_prepare_integration_tl(c_key_self, c_key_incr) &
           bind(c,name='mpas_model_prepare_integration_tl_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_incr  !< Model fields

type(mpas_model), pointer :: self
type(mpas_fields), pointer :: inc

call mpas_model_registry%get(c_key_self, self)
call mpas_fields_registry%get(c_key_incr, inc)

call model_prepare_integration_tl(self, inc)

end subroutine c_mpas_model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine c_mpas_model_propagate(c_key_self, c_key_state) bind(c,name='mpas_model_propagate_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model fields

type(mpas_model), pointer :: self
type(mpas_fields), pointer :: state

call mpas_model_registry%get(c_key_self, self)
call mpas_fields_registry%get(c_key_state,state)

call model_propagate(self, state)

end subroutine c_mpas_model_propagate

! ------------------------------------------------------------------------------

subroutine c_mpas_model_propagate_ad(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='mpas_model_propagate_ad_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_trajectories
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model fields
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure

type(mpas_model),      pointer :: self
type(mpas_fields),  pointer :: inc
type(mpas_trajectory), pointer :: traj

call mpas_model_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_incr,inc)
call mpas_traj_registry%get(c_key_traj,traj)

call model_propagate_ad(self, inc, traj)

end subroutine c_mpas_model_propagate_ad

! ------------------------------------------------------------------------------

subroutine c_mpas_model_propagate_tl(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='mpas_model_propagate_tl_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_trajectories
use mpas_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model fields
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure

type(mpas_model),      pointer :: self
type(mpas_fields),  pointer :: inc
type(mpas_trajectory), pointer :: traj

call mpas_model_registry%get(c_key_self, self)
call mpas_fields_registry%get(c_key_incr,inc)
call mpas_traj_registry%get(c_key_traj,traj)

call model_propagate_tl(self, inc, traj)

end subroutine c_mpas_model_propagate_tl

! ------------------------------------------------------------------------------

subroutine c_mpas_model_prop_traj(c_key_self, c_key_state, c_key_traj) bind(c,name='mpas_model_prop_traj_f90')

use iso_c_binding
use mpas_fields_mod
use mpas_model_mod
use mpas_trajectories

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model fields
integer(c_int), intent(inout) :: c_key_traj  !< Trajectory structure

type(mpas_model),      pointer :: self
type(mpas_fields),      pointer :: state
type(mpas_trajectory), pointer :: traj

call mpas_model_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_state,state)

call mpas_traj_registry%init()            
call mpas_traj_registry%add(c_key_traj)
call mpas_traj_registry%get(c_key_traj,traj)

call model_prop_traj(self, state, traj)

end subroutine c_mpas_model_prop_traj

! ------------------------------------------------------------------------------

subroutine c_mpas_model_wipe_traj(c_key_traj) bind(c,name='mpas_model_wipe_traj_f90')

use iso_c_binding
use mpas_model_mod
use mpas_trajectories

implicit none
integer(c_int), intent(inout)   :: c_key_traj  !< Trajectory structure
type(mpas_trajectory), pointer :: traj

call mpas_traj_registry%get(c_key_traj,traj)
call mpas_traj_registry%remove(c_key_traj)

call model_wipe_traj(traj)

end subroutine c_mpas_model_wipe_traj

! ------------------------------------------------------------------------------
