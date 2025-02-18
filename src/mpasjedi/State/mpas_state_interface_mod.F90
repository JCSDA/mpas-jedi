! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module mpas_state_interface_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use atlas_module, only: atlas_fieldset

!oops
use datetime_mod
use oops_variables_mod

!ufo
use ufo_vars_mod, only: ufo_vars_getindex

!MPAS-Model
use mpas_kind_types, only: StrKIND, RKIND
use mpas_pool_routines, only: mpas_pool_get_config

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_state_mod
use mpas_fields_mod
use mpas_kinds, only : c_real_type

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

type(mpas_fields), pointer :: self
type(mpas_geom), pointer :: geom
type(oops_variables) :: state_vars
type(oops_variables) :: inc_vars
character(len=StrKIND),  pointer :: config_microp_scheme, &
                                    config_radt_cld_scheme
logical,  pointer :: config_microp_re
integer :: ivar, jvar

call mpas_fields_registry%init()
call mpas_fields_registry%add(c_key_self)
call mpas_fields_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)

state_vars = oops_variables(c_state_vars)
inc_vars = oops_variables(c_inc_vars)

call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_radt_cld_scheme', config_radt_cld_scheme)

! TODO(JJG): mpas_re_fields should be shifted to the state variables in the yaml file
if (config_microp_re) then
   do ivar = 1, state_vars % nvars()
      ! Only need re when hydrometeors are in state
      if ( ufo_vars_getindex( mpas_hydrometeor_fields, &
                              state_vars%variable(ivar) ) > 0 ) then
         do jvar = 1, size(mpas_re_fields, 1)
            if (.not. state_vars%has(mpas_re_fields(jvar))) &
               call state_vars%push_back(mpas_re_fields(jvar))
         end do
         exit
      end if
   end do
end if

! TODO(JJG): "nr" should be shifted to the state variables in the yaml file
if (trim(config_microp_scheme) == 'mp_thompson') then
   do ivar = 1, state_vars % nvars()
      ! Only need nr when qr is in state
      if ( trim(state_vars%variable(ivar)) == 'rain_water' ) then
         call state_vars%push_back("rain_number_concentration")
         exit
      end if
   end do
end if

call self%create(geom, state_vars, inc_vars)

end subroutine mpas_state_create_c

! ------------------------------------------------------------------------------

subroutine mpas_state_delete_c(c_key_self) &
      bind(c,name='mpas_state_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)

call self%delete()

call mpas_fields_registry%remove(c_key_self)

end subroutine mpas_state_delete_c

! ------------------------------------------------------------------------------

subroutine mpas_state_zero_c(c_key_self) &
      bind(c,name='mpas_state_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)
call self%zeros()

end subroutine mpas_state_zero_c

! ------------------------------------------------------------------------------

subroutine mpas_state_copy_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_copy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call self%copy(rhs)

end subroutine mpas_state_copy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_axpy_c(c_key_self,c_zz,c_key_rhs) &
      bind(c,name='mpas_state_axpy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_real_type), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
real(kind=RKIND) :: zz

call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)
zz = real(c_zz, kind=RKIND)

call self%axpy(zz,rhs)

end subroutine mpas_state_axpy_c

! ------------------------------------------------------------------------------

subroutine mpas_state_add_incr_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_state_add_incr_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs

call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine mpas_state_add_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_state_read_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_read_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state !< State
type(c_ptr), value, intent(in) :: c_conf  !< Configuration
type(c_ptr), value, intent(in) :: c_dt    !< DateTime

type(mpas_fields), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_state,self)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call self%read_file(f_conf, fdate)

end subroutine mpas_state_read_file_c

! ------------------------------------------------------------------------------

subroutine mpas_state_analytic_init_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_analytic_init_f90')
implicit none
integer(c_int), intent(in) :: c_key_state !< State
type(c_ptr), value, intent(in) :: c_conf  !< Configuration
type(c_ptr), value, intent(in) :: c_dt    !< DateTime

type(mpas_fields), pointer :: state
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_state,state)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call analytic_IC(state, f_conf, fdate)

end subroutine mpas_state_analytic_init_c

! ------------------------------------------------------------------------------

subroutine mpas_state_write_file_c(c_key_state, c_conf, c_dt) &
      bind(c,name='mpas_state_write_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_state !< State
type(c_ptr), value, intent(in) :: c_conf  !< Configuration
type(c_ptr), value, intent(in) :: c_dt    !< DateTime

type(mpas_fields), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_state,self)
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
real(c_real_type), intent(inout) :: pstat(3*kf)

type(mpas_fields), pointer :: self
real(kind=RKIND) :: zstat(3, kf)
integer :: jj, js, jf

call mpas_fields_registry%get(c_key_state,self)

call self%gpnorm(kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = real(zstat(js,jf), kind=c_real_type)
  enddo
enddo

end subroutine mpas_state_gpnorm_c

! ------------------------------------------------------------------------------

subroutine mpas_state_rms_c(c_key_state, prms) &
      bind(c,name='mpas_state_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_state
real(c_real_type), intent(inout) :: prms

type(mpas_fields), pointer :: self
real(kind=RKIND) :: zz

call mpas_fields_registry%get(c_key_state,self)

call self%rms(zz)

prms = real(zz, kind=c_real_type)

end subroutine mpas_state_rms_c

! ------------------------------------------------------------------------------

subroutine mpas_state_sizes_c(c_key_self,nc,nf) &
      bind(c,name='mpas_state_sizes_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nc,nf
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)

nf = self%nf_ci
nc = self%geom%nCellsGlobal

end subroutine mpas_state_sizes_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_state_serial_size_c(c_key_self,c_vsize) &
      bind(c,name='mpas_state_serial_size_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< State
integer(c_size_t),intent(out) :: c_vsize !< Size

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)
call self%serial_size(c_vsize)

end subroutine mpas_state_serial_size_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_state_serialize_c(c_key_self,c_vsize,c_vect_inc) &
      bind(c,name='mpas_state_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self           !< State
integer(c_size_t),intent(in) :: c_vsize           !< Size
real(c_real_type),intent(out) :: c_vect_inc(c_vsize) !< Vector

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)
call self%serialize(c_vsize, c_vect_inc)

end subroutine mpas_state_serialize_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_state_deserialize_c(c_key_self,c_vsize,c_vect_inc,c_index) &
      bind(c,name='mpas_state_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< State
integer(c_size_t),intent(in) :: c_vsize          !< Size
real(c_real_type),intent(in) :: c_vect_inc(c_vsize) !< Vector
integer(c_size_t), intent(inout):: c_index       !< Index

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)
call self%deserialize(c_vsize, c_vect_inc, c_index)

end subroutine mpas_state_deserialize_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_state_to_fieldset_c(c_key_self, c_key_geom, c_vars, c_afieldset, c_include_halo, c_flip_vert_lev) &
 & bind (c,name='mpas_state_to_fieldset_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset
logical(c_bool), intent(in)    :: c_include_halo
logical(c_bool), intent(in)    :: c_flip_vert_lev

type(mpas_fields), pointer :: self
type(mpas_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset
logical :: include_halo
logical :: flip_vert_lev

call mpas_fields_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

include_halo = c_include_halo
flip_vert_lev = c_flip_vert_lev

call self%to_fieldset(geom, vars, afieldset, include_halo, flip_vert_lev)

end subroutine mpas_state_to_fieldset_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_state_from_fieldset_c(c_key_self, c_key_geom, c_vars, c_afieldset, c_include_halo, c_flip_vert_lev) &
 & bind (c,name='mpas_state_from_fieldset_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset
logical(c_bool), intent(in)    :: c_include_halo
logical(c_bool), intent(in)    :: c_flip_vert_lev

type(mpas_fields), pointer :: self
type(mpas_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset
logical :: include_halo
logical :: flip_vert_lev

call mpas_fields_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

include_halo = c_include_halo
flip_vert_lev = c_flip_vert_lev

call self%from_fieldset(geom, vars, afieldset, include_halo, flip_vert_lev)

end subroutine mpas_state_from_fieldset_c

! --------------------------------------------------------------------------------------------------

end module mpas_state_interface_mod

