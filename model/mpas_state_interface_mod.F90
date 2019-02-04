! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module mpas_state_interface_mod

use iso_c_binding
use mpas_state_mod
use mpas_state_utils_mod
use mpas_increment_utils_mod

use mpas_geom_mod
use ufo_vars_mod

use kinds, only: kind_real

!State read/write/init
use datetime_mod

!GetValues+traj
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use mpas_getvaltraj_mod, only: mpas_getvaltraj, mpas_getvaltraj_registry

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!subroutine mpas_state_create_c(c_key_self, c_key_geom, c_key_vars) 
!      bind(c,name='mpas_state_create_f90')
subroutine mpas_state_create_c(c_key_self, c_key_geom, c_vars) &
      bind(c,name='mpas_state_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
!integer(c_int), intent(in) :: c_key_vars !< List of variables
type(c_ptr), intent(in) :: c_vars !< List of variables

type(mpas_state), pointer :: self
type(mpas_geom),  pointer :: geom
type(ufo_vars) :: vars

call mpas_geom_registry%get(c_key_geom, geom)
call mpas_state_registry%init()
call mpas_state_registry%add(c_key_self)
call mpas_state_registry%get(c_key_self,self)

call ufo_vars_setup(vars, c_vars)

call self%create(geom, vars)

end subroutine mpas_state_create_c

! ------------------------------------------------------------------------------

subroutine mpas_state_delete_c(c_key_self) &
      bind(c,name='mpas_state_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_state), pointer :: self

call mpas_state_registry%get(c_key_self,self)

call self%delete()

call mpas_state_registry%remove(c_key_self)

end subroutine mpas_state_delete_c

! ------------------------------------------------------------------------------

subroutine mpas_state_zero_c(c_key_self) &
      bind(c,name='mpas_state_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_state), pointer :: self

call mpas_state_registry%get(c_key_self,self)
call self%zeros()

end subroutine mpas_state_zero_c

! ------------------------------------------------------------------------------

subroutine mpas_state_copy_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_copy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_state), pointer :: self
type(mpas_state), pointer :: rhs
call mpas_state_registry%get(c_key_self,self)
call mpas_state_registry%get(c_key_rhs,rhs)

call self%copy(rhs)

end subroutine mpas_state_copy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_axpy_c(c_key_self,c_zz,c_key_rhs) &
      bind(c,name='mpas_state_axpy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_state), pointer :: self
type(mpas_state), pointer :: rhs
real(kind=kind_real) :: zz

call mpas_state_registry%get(c_key_self,self)
call mpas_state_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%axpy(zz,rhs)

end subroutine mpas_state_axpy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_add_incr_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_add_incr_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(mpas_state), pointer :: self
type(mpas_increment), pointer :: rhs

call mpas_state_registry%get(c_key_self,self)
call mpas_increment_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine mpas_state_add_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_state_change_resol_c(c_key_state,c_key_rhs) &
      bind(c,name='mpas_state_change_resol_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: c_key_rhs
type(mpas_state), pointer :: state, rhs

call mpas_state_registry%get(c_key_state,state)
call mpas_state_registry%get(c_key_rhs,rhs)

call state%change_resol(rhs)

end subroutine mpas_state_change_resol_c

! ------------------------------------------------------------------------------

subroutine mpas_state_read_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_read_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(mpas_state), pointer :: self
type(datetime) :: fdate

call mpas_state_registry%get(c_key_state,self)
call c_f_datetime(c_dt, fdate)
call self%read_file(c_conf, fdate)

end subroutine mpas_state_read_file_c

! ------------------------------------------------------------------------------

subroutine mpas_state_analytic_init_c(c_key_state, c_key_geom, c_conf, c_dt) &
      bind(c,name='mpas_state_analytic_init_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(mpas_state), pointer :: state
type(mpas_geom), pointer :: geom
type(datetime) :: fdate

call mpas_state_registry%get(c_key_state,state)
call mpas_geom_registry%get(c_key_geom,geom)
call c_f_datetime(c_dt, fdate)
call analytic_IC(state, geom, c_conf, fdate)

end subroutine mpas_state_analytic_init_c

! ------------------------------------------------------------------------------

subroutine mpas_state_write_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_write_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime

type(mpas_state), pointer :: self
type(datetime) :: fdate

call mpas_state_registry%get(c_key_state,self)
call c_f_datetime(c_dt, fdate)
call self%write_file( c_conf, fdate)

end subroutine mpas_state_write_file_c

! ------------------------------------------------------------------------------

subroutine mpas_state_gpnorm_c(c_key_state, kf, pstat) &
      bind(c,name='mpas_state_gpnorm_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(mpas_state), pointer :: self
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call mpas_state_registry%get(c_key_state,self)

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

type(mpas_state), pointer :: self
real(kind=kind_real) :: zz

call mpas_state_registry%get(c_key_state,self)

call self%rms(zz)

prms = zz

end subroutine mpas_state_rms_c

! ------------------------------------------------------------------------------

subroutine mpas_state_getvalues_notraj_c(c_key_state,c_key_loc,c_vars,c_key_gom) &
      bind(c,name='mpas_state_getvalues_notraj_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in) :: c_vars  !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
type(mpas_state), pointer :: state
type(ufo_locs),  pointer :: locs
type(ufo_vars)  :: vars
type(ufo_geovals),  pointer :: gom

call mpas_state_registry%get(c_key_state, state)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_vars_setup(vars, c_vars)
call ufo_geovals_registry%get(c_key_gom, gom)

call getvalues(state, locs, vars, gom)

end subroutine mpas_state_getvalues_notraj_c

! ------------------------------------------------------------------------------

subroutine mpas_state_getvalues_c(c_key_state,c_key_loc,c_vars,c_key_gom,c_key_traj) &
      bind(c,name='mpas_state_getvalues_f90')
implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in) :: c_vars  !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in), optional :: c_key_traj !< Trajectory for interpolation/transforms
type(mpas_state), pointer :: state
type(ufo_locs),  pointer :: locs
type(ufo_vars)  :: vars
type(ufo_geovals),  pointer :: gom
type(mpas_getvaltraj), pointer :: traj

call mpas_state_registry%get(c_key_state, state)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_vars_setup(vars, c_vars)
call ufo_geovals_registry%get(c_key_gom, gom)
call mpas_getvaltraj_registry%get(c_key_traj, traj)

call getvalues(state, locs, vars, gom, traj)

end subroutine mpas_state_getvalues_c

! ------------------------------------------------------------------------------

subroutine mpas_state_sizes_c(c_key_self,nc,nf) &
      bind(c,name='mpas_state_sizes_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nc,nf
type(mpas_state), pointer :: self

call mpas_state_registry%get(c_key_self,self)

nf = self%nf_ci
nc = self%geom%nCellsGlobal

end subroutine mpas_state_sizes_c

! ------------------------------------------------------------------------------   

end module mpas_state_interface_mod

