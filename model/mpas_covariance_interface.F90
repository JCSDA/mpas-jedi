! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_mpas_b_setup(c_key_self, c_conf, c_key_geom) &
          & bind (c,name='mpas_b_setup_f90')

use iso_c_binding
use mpas_covariance_mod
use mpas_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure
type(c_ptr), intent(in)    :: c_conf         !< Configuration
integer(c_int), intent(in) :: c_key_geom     !< Geometry
type(mpas_covar), pointer :: self
type(mpas_geom),  pointer :: geom

call mpas_geom_registry%get(c_key_geom, geom)
call mpas_covar_registry%init()
call mpas_covar_registry%add(c_key_self)
call mpas_covar_registry%get(c_key_self, self)

call mpas_covar_setup(self, geom, c_conf)

end subroutine c_mpas_b_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_b_delete(c_key_self) bind (c,name='mpas_b_delete_f90')

use iso_c_binding
use mpas_covariance_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure
type(mpas_covar), pointer :: self

call mpas_covar_registry%get(c_key_self,self)
call mpas_covar_delete(self)
call mpas_covar_registry%remove(c_key_self)

end subroutine c_mpas_b_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse of covariance

subroutine c_mpas_b_inv_mult(c_key_self, c_key_in, c_key_out) bind(c,name='mpas_b_invmult_f90')

use iso_c_binding
use mpas_covariance_mod
use mpas_fields_mod
use kinds
use mpas_framework !BJJ

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out
type(mpas_covar), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_covar_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)
!TODO BJJ
!Implement this
!xout = xin
!xout => xin
   call mpas_pool_empty_pool(xout % subFields)
   call mpas_pool_destroy_pool(xout % subFields)
   xout % nf = xin % nf
   call mpas_pool_create_pool(xout % subFields, xout % nf)
   call mpas_pool_clone_pool(xin % subFields, xout % subFields)


!TODO BJJ
!Implement this
!xout = xin
!xout => xin
   call mpas_pool_empty_pool(xout % subFields)
   call mpas_pool_destroy_pool(xout % subFields)
   xout % nf = xin % nf
   call mpas_pool_create_pool(xout % subFields, xout % nf)
   call mpas_pool_clone_pool(xin % subFields, xout % subFields)

!call mpas_covar_sqrt_inv_mult(self%nx,self%ny,xctl,xin,self)
!call zeros(xout)
!call mpas_covar_sqrt_inv_mult_ad(self%nx,self%ny,xctl,xout,self)

end subroutine c_mpas_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by covariance

subroutine c_mpas_b_mult(c_key_self, c_key_in, c_key_out) bind(c,name='mpas_b_mult_f90')

use iso_c_binding
use mpas_covariance_mod
use mpas_fields_mod
use kinds
use mpas_framework !BJJ

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out
type(mpas_covar), pointer :: self
type(mpas_field), pointer :: xin
type(mpas_field), pointer :: xout

call mpas_covar_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_in,xin)
call mpas_field_registry%get(c_key_out,xout)

   write(*,*) '---- inside sub c_mpas_b_mult ----'
!TODO BJJ
!Implement this
!xout = xin
!xout => xin
   call mpas_pool_empty_pool(xout % subFields)
   call mpas_pool_destroy_pool(xout % subFields)
   xout % nf = xin % nf
   call mpas_pool_create_pool(xout % subFields, xout % nf)
   call mpas_pool_clone_pool(xin % subFields, xout % subFields)

!call mpas_covar_sqrt_mult_ad(self%nx,self%ny,xin,xctl,self)
!call zeros(xout)
!call mpas_covar_sqrt_mult(self%nx,self%ny,xout,xctl,self)

end subroutine c_mpas_b_mult

! ------------------------------------------------------------------------------

!> Generate randomized increment

subroutine c_mpas_b_randomize(c_key_self, c_key_out) bind(c,name='mpas_b_randomize_f90')

use iso_c_binding
use mpas_covariance_mod
use mpas_fields_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_out
type(mpas_covar), pointer :: self
type(mpas_field), pointer :: xout

call mpas_covar_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_out,xout)

call random(xout)

end subroutine c_mpas_b_randomize

! ------------------------------------------------------------------------------
