! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module mpas_increment_interface_mod

use atlas_module, only: atlas_fieldset
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod

!mpas-jedi
use mpas_geom_mod
use mpas_increment_mod
use mpas_fields_mod

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine mpas_increment_create_c(c_key_self, c_key_geom, c_vars) &
      bind(c,name='mpas_increment_create_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(c_ptr), value, intent(in) :: c_vars !< List of variables

type(mpas_fields), pointer :: self
type(mpas_geom), pointer :: geom
type(oops_variables) :: vars

call mpas_geom_registry%get(c_key_geom, geom)
call mpas_fields_registry%init()
call mpas_fields_registry%add(c_key_self)
call mpas_fields_registry%get(c_key_self,self)

vars = oops_variables(c_vars)
call self%create(geom, vars, vars)

end subroutine mpas_increment_create_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_delete_c(c_key_self) &
      bind(c,name='mpas_increment_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)

call self%delete()

call mpas_fields_registry%remove(c_key_self)

end subroutine mpas_increment_delete_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_zero_c(c_key_self) &
      bind(c,name='mpas_increment_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)
call self%zeros()

end subroutine mpas_increment_zero_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_ones_c(c_key_self) &
      bind(c,name='mpas_increment_ones_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)
call self%ones()

end subroutine mpas_increment_ones_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_dirac_c(c_key_self,c_conf) &
      bind(c,name='mpas_increment_dirac_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_conf !< Configuration

type(mpas_fields), pointer :: self
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_self,self)
f_conf = fckit_configuration(c_conf)
call dirac(self,f_conf)

end subroutine mpas_increment_dirac_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_random_c(c_key_self) &
      bind(c,name='mpas_increment_random_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)
call self%random()

end subroutine mpas_increment_random_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_copy_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_increment_copy_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call self%copy(rhs)

end subroutine mpas_increment_copy_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_self_add_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_increment_self_add_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call self%self_add(rhs)

end subroutine mpas_increment_self_add_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_self_schur_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_increment_self_schur_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call self%self_schur(rhs)

end subroutine mpas_increment_self_schur_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_self_sub_c(c_key_self,c_key_rhs) &
      bind(c,name='mpas_increment_self_sub_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)

call self%self_sub(rhs)

end subroutine mpas_increment_self_sub_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_self_mul_c(c_key_self,c_zz) &
      bind(c,name='mpas_increment_self_mul_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(mpas_fields), pointer :: self
real(kind=kind_real) :: zz

call mpas_fields_registry%get(c_key_self,self)
zz = c_zz

call self%self_mult(zz)

end subroutine mpas_increment_self_mul_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_axpy_inc_c(c_key_self,c_zz,c_key_rhs) &
      bind(c,name='mpas_increment_axpy_inc_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
real(kind=kind_real) :: zz

call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%axpy(zz,rhs)

end subroutine mpas_increment_axpy_inc_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_axpy_state_c(c_key_self,c_zz,c_key_rhs) &
      bind(c,name='mpas_increment_axpy_state_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_fields), pointer :: self
type(mpas_fields), pointer :: rhs
real(kind=kind_real) :: zz

call mpas_fields_registry%get(c_key_self,self)
call mpas_fields_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%axpy(zz,rhs)

end subroutine mpas_increment_axpy_state_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_dot_prod_c(c_key_inc1,c_key_inc2,c_prod) &
      bind(c,name='mpas_increment_dot_prod_f90')
implicit none
integer(c_int), intent(in)    :: c_key_inc1, c_key_inc2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(mpas_fields), pointer :: inc1, inc2

call mpas_fields_registry%get(c_key_inc1,inc1)
call mpas_fields_registry%get(c_key_inc2,inc2)

call inc1%dot_prod(inc2,zz)

c_prod = zz

end subroutine mpas_increment_dot_prod_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) &
      bind(c,name='mpas_increment_diff_incr_f90')
implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
type(mpas_fields), pointer :: lhs
type(mpas_fields), pointer :: x1
type(mpas_fields), pointer :: x2

call mpas_fields_registry%get(c_key_lhs,lhs)
call mpas_fields_registry%get(c_key_x1,x1)
call mpas_fields_registry%get(c_key_x2,x2)

call diff_incr(lhs,x1,x2)

end subroutine mpas_increment_diff_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_change_resol_c(c_key_inc,c_key_rhs) &
      bind(c,name='mpas_increment_change_resol_f90')
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_rhs
type(mpas_fields), pointer :: inc, rhs

call mpas_fields_registry%get(c_key_inc,inc)
call mpas_fields_registry%get(c_key_rhs,rhs)

call inc%change_resol(rhs)

end subroutine mpas_increment_change_resol_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_set_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
      bind (c,name='mpas_increment_set_atlas_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset
type(mpas_fields), pointer :: self
type(mpas_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call mpas_fields_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call set_atlas(self, geom, vars, afieldset)

end subroutine mpas_increment_set_atlas_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_to_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
      bind (c,name='mpas_increment_to_atlas_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset
type(mpas_fields), pointer :: self
type(mpas_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call mpas_fields_registry%get(c_key_self,self)
call mpas_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call to_atlas(self, geom, vars, afieldset)

end subroutine mpas_increment_to_atlas_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_from_atlas_c(c_key_self,c_key_geom,c_vars,c_afieldset) &
      bind (c,name='mpas_increment_from_atlas_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset
type(mpas_fields), pointer :: self
type(mpas_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call mpas_fields_registry%get(c_key_self, self)
call mpas_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call from_atlas(self, geom, vars, afieldset)

end subroutine mpas_increment_from_atlas_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_read_file_c(c_key_inc, c_conf, c_dt) &
      bind(c,name='mpas_increment_read_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_inc  !< Fields
type(c_ptr), value, intent(in) :: c_conf !< Configuration
type(c_ptr), value, intent(in) :: c_dt   !< DateTime

type(mpas_fields), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_inc,self)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call self%read_file(f_conf, fdate)

end subroutine mpas_increment_read_file_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_write_file_c(c_key_inc, c_conf, c_dt) &
      bind(c,name='mpas_increment_write_file_f90')
implicit none
integer(c_int), intent(in) :: c_key_inc  !< Fields
type(c_ptr), value, intent(in) :: c_conf !< Configuration
type(c_ptr), value, intent(in) :: c_dt   !< DateTime

type(mpas_fields), pointer :: self
type(datetime) :: fdate
type(fckit_configuration) :: f_conf

call mpas_fields_registry%get(c_key_inc,self)
call c_f_datetime(c_dt, fdate)
f_conf = fckit_configuration(c_conf)
call self%write_file(f_conf, fdate)

end subroutine mpas_increment_write_file_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_gpnorm_c(c_key_inc, kf, pstat) &
      bind(c,name='mpas_increment_gpnorm_f90')
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(mpas_fields), pointer :: self
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call mpas_fields_registry%get(c_key_inc,self)

call self%gpnorm(kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine mpas_increment_gpnorm_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_rms_c(c_key_inc, prms) &
      bind(c,name='mpas_increment_rms_f90')
implicit none
integer(c_int), intent(in) :: c_key_inc
real(c_double), intent(inout) :: prms

type(mpas_fields), pointer :: self
real(kind=kind_real) :: zz

call mpas_fields_registry%get(c_key_inc,self)

call self%rms(zz)

prms = zz

end subroutine mpas_increment_rms_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_print_c(c_key_self) &
      bind(c,name='mpas_increment_print_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)

!call increment_print(self)

end subroutine mpas_increment_print_c

! ------------------------------------------------------------------------------

subroutine mpas_increment_sizes_c(c_key_self,nc,nf) &
      bind(c,name='mpas_increment_sizes_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nc,nf
type(mpas_fields), pointer :: self

call mpas_fields_registry%get(c_key_self,self)

nf = self%nf_ci
nc = self%geom%nCellsGlobal

end subroutine mpas_increment_sizes_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_increment_serial_size_c(c_key_self,c_vsize) &
      bind(c,name='mpas_increment_serial_size_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Increment
integer(c_size_t),intent(out) :: c_vsize !< Size

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)
call self%serial_size(c_vsize)

end subroutine mpas_increment_serial_size_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_increment_serialize_c(c_key_self,c_vsize,c_vect_inc) &
      bind(c,name='mpas_increment_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self           !< Increment
integer(c_size_t),intent(in) :: c_vsize           !< Size
real(c_double),intent(out) :: c_vect_inc(c_vsize) !< Vector

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)
! Call Fortran
call self%serialize(c_vsize, c_vect_inc)

end subroutine mpas_increment_serialize_c

! --------------------------------------------------------------------------------------------------

subroutine mpas_increment_deserialize_c(c_key_self,c_vsize,c_vect_inc,c_index) &
      bind(c,name='mpas_increment_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Increment
integer(c_size_t),intent(in) :: c_vsize          !< Size
real(c_double),intent(in) :: c_vect_inc(c_vsize) !< Vector
integer(c_size_t), intent(inout):: c_index       !< Index

type(mpas_fields),pointer :: self

call mpas_fields_registry%get(c_key_self, self)

! Call Fortran
call self%deserialize(c_vsize, c_vect_inc, c_index)

end subroutine mpas_increment_deserialize_c
! ------------------------------------------------------------------------------

end module mpas_increment_interface_mod

