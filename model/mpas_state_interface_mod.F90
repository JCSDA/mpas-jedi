! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module mpas_state_interface_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real
use oops_variables_mod

!mpas-jedi
use mpas_geom_mod
use mpas_state_mod
use mpas_field_utils_mod
use mpas_getvaltraj_mod, only: mpas_getvaltraj, mpas_getvaltraj_registry

!MPAS-Model
use mpas_kind_types, only: StrKIND
use mpas_pool_routines, only: mpas_pool_get_config

!State read/write/init
use datetime_mod

!GetValues+traj
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine mpas_state_create_c(c_key_self, c_key_geom, c_state_vars, c_inc_vars) &
      bind(c,name='mpas_state_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom       !< Geometry
type(c_ptr), value, intent(in) :: c_state_vars !< List of state variables
type(c_ptr), value, intent(in) :: c_inc_vars   !< List of increment variables

type(mpas_field), pointer :: self
type(mpas_geom), pointer :: geom
type(oops_variables) :: state_vars
type(oops_variables) :: inc_vars
character(len=StrKIND),  pointer :: config_microp_scheme
logical,  pointer :: config_microp_re

call mpas_field_registry%init()
call mpas_field_registry%add(c_key_self)
call mpas_field_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)

state_vars = oops_variables(c_state_vars)
inc_vars = oops_variables(c_inc_vars)

call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
write(*,*) "config_microp_re=",config_microp_re
write(*,*) "config_microp_scheme=",config_microp_scheme
if (config_microp_re) then
   call state_vars%push_back(mpas_re_fields)
end if
if (trim(config_microp_scheme) == 'mp_thompson') then
   call state_vars%push_back("index_nr")
end if
call self%create(geom, state_vars, inc_vars)

end subroutine mpas_state_create_c

! ------------------------------------------------------------------------------

subroutine mpas_state_delete_c(c_key_self) &
      bind(c,name='mpas_state_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)

call self%delete()

call mpas_field_registry%remove(c_key_self)

end subroutine mpas_state_delete_c

! ------------------------------------------------------------------------------

subroutine mpas_state_zero_c(c_key_self) &
      bind(c,name='mpas_state_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)
call self%zeros()

end subroutine mpas_state_zero_c

! ------------------------------------------------------------------------------

subroutine mpas_state_copy_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_copy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call self%copy(rhs)

end subroutine mpas_state_copy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_axpy_c(c_key_self,c_zz,c_key_rhs) &
      bind(c,name='mpas_state_axpy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
real(kind=kind_real) :: zz

call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%axpy(zz,rhs)

end subroutine mpas_state_axpy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_add_incr_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_add_incr_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs

call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine mpas_state_add_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_state_change_resol_c(c_key_state,c_key_rhs) &
      bind(c,name='mpas_state_change_resol_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: c_key_rhs
type(mpas_field), pointer :: state, rhs

call mpas_field_registry%get(c_key_state,state)
call mpas_field_registry%get(c_key_rhs,rhs)

call state%change_resol(rhs)

end subroutine mpas_state_change_resol_c

! ------------------------------------------------------------------------------

subroutine mpas_state_read_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_read_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(mpas_field), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_field_registry%get(c_key_state,self)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call self%read_file(f_conf, fdate)

end subroutine mpas_state_read_file_c

! ------------------------------------------------------------------------------

subroutine mpas_state_analytic_init_c(c_key_state, c_key_geom, c_conf, c_dt) &
      bind(c,name='mpas_state_analytic_init_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(mpas_field), pointer :: state
type(mpas_geom), pointer :: geom
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_field_registry%get(c_key_state,state)
call mpas_geom_registry%get(c_key_geom,geom)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call analytic_IC(state, geom, f_conf, fdate)

end subroutine mpas_state_analytic_init_c

! ------------------------------------------------------------------------------

subroutine mpas_state_write_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_write_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime

type(mpas_field), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_field_registry%get(c_key_state,self)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call self%write_file( f_conf, fdate)

end subroutine mpas_state_write_file_c

! ------------------------------------------------------------------------------

subroutine mpas_state_gpnorm_c(c_key_state, kf, pstat) &
      bind(c,name='mpas_state_gpnorm_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(mpas_field), pointer :: self
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call mpas_field_registry%get(c_key_state,self)

call self%gpnorm(kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine mpas_state_gpnorm_c

! ------------------------------------------------------------------------------

subroutine mpas_state_rms_c(c_key_state, prms) &
      bind(c,name='mpas_state_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
real(c_double), intent(inout) :: prms

type(mpas_field), pointer :: self
real(kind=kind_real) :: zz

call mpas_field_registry%get(c_key_state,self)

call self%rms(zz)

prms = zz

end subroutine mpas_state_rms_c

! ------------------------------------------------------------------------------

subroutine mpas_state_getvalues_notraj_c(c_key_state,c_key_loc,c_vars,c_key_gom) &
      bind(c,name='mpas_state_getvalues_notraj_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), value, intent(in) :: c_vars  !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values

type(mpas_field),   pointer :: state
type(ufo_locs),     pointer :: locs
type(oops_variables)        :: vars
type(ufo_geovals),  pointer :: gom

call mpas_field_registry%get(c_key_state, state)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)

vars = oops_variables(c_vars)
call getvalues(state, locs, vars, gom)

end subroutine mpas_state_getvalues_notraj_c

! ------------------------------------------------------------------------------

subroutine mpas_state_getvalues_c(c_key_state,c_key_loc,c_vars,c_key_gom,c_key_traj) &
      bind(c,name='mpas_state_getvalues_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), value, intent(in) :: c_vars  !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in), optional :: c_key_traj !< Trajectory for interpolation/transforms

type(mpas_field), pointer :: state
type(ufo_locs), pointer :: locs
type(oops_variables) :: vars
type(ufo_geovals), pointer :: gom
type(mpas_getvaltraj), pointer :: traj

call mpas_field_registry%get(c_key_state, state)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call mpas_getvaltraj_registry%get(c_key_traj, traj)

vars = oops_variables(c_vars)
call getvalues(state, locs, vars, gom, traj)

end subroutine mpas_state_getvalues_c

! ------------------------------------------------------------------------------

subroutine mpas_state_sizes_c(c_key_self,nc,nf) &
      bind(c,name='mpas_state_sizes_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nc,nf
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)

nf = self%nf_ci
nc = self%geom%nCellsGlobal

end subroutine mpas_state_sizes_c

! ------------------------------------------------------------------------------   

end module mpas_state_interface_mod

